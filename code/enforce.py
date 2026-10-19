"""Process-isolation runtime that builds the sole-egress binding.

Motivation: `governance.EgressMonitor` is a plaintext library call, so a study
author who does not call `release` writes withheld quantities straight out. This
module runs the study computation as a CHILD process whose ONLY writable channel
is one pre-opened pipe descriptor wired to the parent's monitor, and installs a
PEP-578 audit hook that denies write-opens, sockets, subprocesses and ctypes.

Threat model = the A2 negligent/careless study-code author (ordinary Python, no
native code, not root), matching the paper. The parent runtime and a local
user with process control stay in the TCB.

Not an MPC implementation; a focused egress-binding prototype.
"""
from __future__ import annotations

import json
import os
import sys
import numpy as np

from common import simulate_population, split_population3, release_objects, \
    allele_freq, af_lrt_attack, auroc
from governance import make_manifest, EgressMonitor, Ledger, Violation, POLICIES

SCRATCH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "results", "scratch")
os.makedirs(SCRATCH, exist_ok=True)


class SandboxViolation(Exception):
    pass


WRITE_BITS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND


def install_egress_audit_hook(allowed_write_fds):
    """Deny every natural egress route an ordinary study script uses, except
    writes through fds in allowed_write_fds (the sanctioned pipe)."""
    allowed = set(int(f) for f in allowed_write_fds)

    def _hook(event, args):
        if event == "open":
            path, mode, flags = args
            if isinstance(mode, str) and any(c in mode for c in "wxa+"):
                raise SandboxViolation(f"open({path!r},{mode!r})")
            if isinstance(flags, int) and (flags & WRITE_BITS):
                raise SandboxViolation(f"open({path!r},flags={flags})")
        elif event in ("socket.connect", "socket.bind", "socket.getaddrinfo",
                       "socket.gethostbyname", "socket.__new__", "socket.sendto",
                       "socket.sendmsg", "socket.send",
                       "subprocess.Popen", "os.system", "os.exec", "os.spawn",
                       "ctypes.dlopen", "ctypes.dlsym"):
            raise SandboxViolation(event)

    sys.addaudithook(_hook)


def _jsonable(rel):
    """The careless author dumps EVERY quantity it has, declared or not."""
    def arr(a):
        return np.asarray(a).tolist()

    out = dict(af=arr(rel["af"]), eig=arr(rel["eig"]),
               lambda_gc=float(rel["lambda_gc"]), loadings=arr(rel["loadings"]),
               missing=arr(rel["missing"]), pvals=arr(rel["pvals"]),
               af_by_group={str(k): arr(v) for k, v in rel["af_by_group"].items()})
    return out


def child_main():
    """Child entrypoint: install the hook, optionally exercise the bypass, run
    the study computation, write one JSON blob to the sanctioned fd, exit."""
    out_fd = int(os.environ["OUT_FD"])
    seed = int(os.environ.get("SEED", "0"))
    inject = bool(int(os.environ.get("INJECT", "0")))
    inject_native = bool(int(os.environ.get("INJECT_NATIVE", "0")))
    leak_path = os.environ.get("LEAK_PATH", "")
    report = dict(seed=seed)

    install_egress_audit_hook(allowed_write_fds=[out_fd])

    if inject and leak_path:
        try:
            with open(leak_path, "wb") as f:            # natural egress route
                f.write(b"loadings")
            report["natural_egress_blocked"] = False
        except SandboxViolation as e:
            report["natural_egress_blocked"] = True
            report["natural_egress_reason"] = str(e)

    if inject_native:
        # the threat model names "no native code": attempt the native-code route
        # a study author could use to write past the audit hook.
        try:
            import ctypes
            libc = ctypes.CDLL(None)                    # audit hook: ctypes.dlopen
            report["native_egress_blocked"] = False
            report["native_egress_reason"] = "ctypes.CDLL(None) returned (not blocked)"
        except SandboxViolation as e:
            report["native_egress_blocked"] = True
            report["native_egress_reason"] = str(e)

    # UDP egress route: ordinary Python can create a datagram socket and sendto
    # without any connect/bind. Test it explicitly (the hook denies the socket
    # creation itself when enforced).
    try:
        import socket as _socket
        _s = _socket.socket(_socket.AF_INET, _socket.SOCK_DGRAM)
        _s.sendto(b"loadings", ("127.0.0.1", 9))
        report["udp_egress_blocked"] = False
        report["udp_egress_reason"] = "UDP sendto returned (not blocked)"
    except SandboxViolation as e:
        report["udp_egress_blocked"] = True
        report["udp_egress_reason"] = str(e)
    except Exception as e:
        report["udp_egress_blocked"] = None
        report["udp_egress_reason"] = f"other:{type(e).__name__}"

    pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
    site, ref, non = split_population3(pop, 300, 300, 300)
    rel = release_objects(site)
    report["emitted"] = _jsonable(rel)
    data = json.dumps(report).encode()
    os.write(out_fd, data)
    os.close(out_fd)


def run_isolated_study(manifest, ledger, seed=0, inject=False, inject_native=False):
    """Parent: spawn the child with a single writable pipe fd wired to the
    monitor; filter every object that arrives; return the outcome."""
    leak_path = None
    if inject:
        leak_path = os.path.join(SCRATCH, f"_s39_iso_leak_{seed}.bin")
        if os.path.exists(leak_path):
            os.remove(leak_path)
    r, w = os.pipe()
    code_dir = os.path.dirname(os.path.abspath(__file__))
    env = dict(os.environ, OUT_FD=str(w), SEED=str(seed),
               INJECT=str(int(inject)), INJECT_NATIVE=str(int(inject_native)),
               LEAK_PATH=leak_path or "",
               PYTHONPATH=code_dir, PYTHONDONTWRITEBYTECODE="1")
    import subprocess
    proc = subprocess.Popen([sys.executable, "-c",
                             "import enforce; enforce.child_main()"],
                            env=env, pass_fds=(w,), cwd=code_dir,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.close(w)
    # read the whole payload first (child may exceed the pipe buffer), then wait
    chunks = []
    with os.fdopen(r, "rb") as f:
        chunks.append(f.read())
    rc = proc.wait()
    raw = b"".join(chunks)

    if rc != 0 or not raw:
        # the child's stdout/stderr are /dev/null (no observable channel besides
        # the sanctioned pipe); a non-zero exit is reported by status alone.
        return dict(seed=seed, child_failed=True, exit=rc,
                    natural_egress_blocked=None, sanctioned_undeclared_blocked=None,
                    declared_ok=False, chain_verified=False,
                    governed_af_auc=float("nan"), governed_pc_auc=float("nan"))

    rep = json.loads(raw.decode())
    blob = rep["emitted"]
    # restore numeric types the monitor expects
    blob = {k: (np.asarray(v) if isinstance(v, list) else v) for k, v in blob.items()}
    blob["af_by_group"] = {int(k): np.asarray(v) for k, v in blob["af_by_group"].items()}

    mon = EgressMonitor(manifest, ledger)
    sanctioned_blocked = False
    try:
        mon.release(blob, "loadings")               # undeclared -> must block
    except Violation:
        sanctioned_blocked = True
    declared = {}
    for qid in ("af", "eig", "lambda_gc"):
        declared[qid] = mon.release(blob, qid)
    chain_ok = ledger.verify_chain()

    # score the transported governed release on a fresh deterministic panel
    pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
    site, ref, non = split_population3(pop, 300, 300, 300)
    rel = release_objects(site)
    mask = rel["mask"]
    af_r = allele_freq(ref["G"])
    y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])
    am = af_lrt_attack(site["G"], declared["af"], af_r, mask=mask)
    an = af_lrt_attack(non["G"], declared["af"], af_r, mask=mask)
    af_auc = float(auroc(np.concatenate([am, an]), y))
    rng = np.random.default_rng(seed)
    Lr = rng.normal(size=rel["loadings"].shape); Lr, _ = np.linalg.qr(Lr.T); Lr = Lr.T
    mean_s = site["G"][:, mask].mean(0)

    def pc_dist(G, load):
        X = G[:, mask].astype(np.float64) - mean_s[None, :]
        return np.linalg.norm((X @ load[:, mask].T) / np.sqrt(max(1, mask.sum())), axis=1)
    dm = pc_dist(site["G"], Lr); dn = pc_dist(non["G"], Lr)
    a = float(auroc(np.concatenate([dm, dn]), y))
    pc_auc = max(a, 1 - a)

    if leak_path and os.path.exists(leak_path):
        os.remove(leak_path)
    return dict(seed=seed, child_failed=False,
                natural_egress_blocked=rep.get("natural_egress_blocked"),
                natural_egress_reason=rep.get("natural_egress_reason", ""),
                native_egress_blocked=rep.get("native_egress_blocked"),
                native_egress_reason=rep.get("native_egress_reason", ""),
                udp_egress_blocked=rep.get("udp_egress_blocked"),
                udp_egress_reason=rep.get("udp_egress_reason", ""),
                sanctioned_undeclared_blocked=bool(sanctioned_blocked),
                declared_ok=bool(np.allclose(declared["af"], np.round(blob["af"], 1))),
                chain_verified=bool(chain_ok),
                n_ledger_records=len(ledger.records),
                governed_af_auc=af_auc, governed_pc_auc=pc_auc)


def bypass_suite(seed=0, native=False):
    """Run the isolated arm WITH the injection; returns the isolation outcome."""
    return run_isolated_study(make_manifest(POLICIES["governed"], f"inj{seed}"),
                              Ledger(), seed=seed, inject=True, inject_native=native)


if __name__ == "__main__":
    import pprint
    pprint.pprint(bypass_suite(0))
