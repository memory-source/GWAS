"""S39 - built sole-egress binding.

The governance prototype's EgressMonitor is a plaintext library call: a study
author who simply does not call `release` writes the withheld quantity straight
out. The paper's Discussion admits the "sole-egress binding" is argued, not
built. This script (A) demonstrates the bypass as a positive control and (B)
closes it with a minimal process-isolation runtime in code/enforce.py:

  - the study computation runs as a CHILD process whose ONLY writable channel is
    a single pre-opened pipe descriptor wired to the parent's EgressMonitor;
  - a startup PEP-578 audit hook denies write-opens, sockets, subprocesses and
    ctypes loading;
  - the parent treats the pipe as untrusted and filters every object through the
    manifest monitor (second enforcement point).

Threat model = the A2 negligent/careless study-code author (ordinary Python,
no native code, not root), matching the paper. Out of scope: a local
user with process control, and the parent runtime (both stay in the TCB).

Output: results/raw/s39_enforcement.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import simulate_population, split_population3, release_objects
from governance import make_manifest, EgressMonitor, Ledger, Violation, POLICIES
import enforce as E
import s23_governed as S23

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def library_arm(seed=0):
    """Positive control A: the SAME exfiltration line against the plaintext
    prototype (no runtime) -- the withheld loadings land on disk."""
    pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
    site, ref, non = split_population3(pop, 300, 300, 300)
    rel = release_objects(site)                       # full-precision object
    led = Ledger(); mon = EgressMonitor(make_manifest(POLICIES["governed"], "lib"), led)
    leak_path = os.path.join(OUT, f"_s39_leak_library_seed{seed}.bin")
    if os.path.exists(leak_path):
        os.remove(leak_path)
    # the study author never calls mon.release(rel, "loadings"); instead:
    with open(leak_path, "wb") as f:
        f.write(rel["loadings"].tobytes())
    wrote = os.path.exists(leak_path)
    if wrote:
        os.remove(leak_path)
    return dict(seed=seed, loadings_exfiltrated=bool(wrote))


def pergroup_delivery(seed=0):
    """A per-group (dict-valued) AF object is a normal governed release object;
    check it flows through the SAME monitor at full precision and at 2dp, with
    the ledger chained, instead of raising."""
    pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
    site, _, _ = split_population3(pop, 300, 300, 300)
    rel = release_objects(site)
    out = {}
    for tag, prec in (("full", "full"), ("2dp", {"kind": "decimals", "d": 2})):
        pol = {q: dict(released=False) for q in ("af", "eig", "lambda_gc")}
        pol["af_by_group"] = dict(released=True, precision=prec)
        mon = EgressMonitor(make_manifest(pol, f"pg{tag}{seed}"), Ledger())
        got = mon.release(rel, "af_by_group")
        ok = isinstance(got, dict) and all(hasattr(v, "__len__") for v in got.values())
        out[tag] = dict(delivered=bool(ok), chained=bool(mon.ledger.verify_chain()),
                        n_groups=len(got) if ok else 0)
    return out


def run(seeds=(0, 1, 2)):
    per_seed = []
    for s in seeds:
        lib = library_arm(s)
        iso = E.bypass_suite(seed=s)
        native = E.bypass_suite(seed=s, native=True)
        transp = E.run_isolated_study(make_manifest(POLICIES["governed"], f"iso{s}"),
                                      Ledger(), seed=s)
        pg = pergroup_delivery(s)
        per_seed.append(dict(seed=s,
                             library_arm=dict(loadings_exfiltrated=lib["loadings_exfiltrated"]),
                             isolated=iso, native=native, transport=transp,
                             pergroup=pg))
        print("seed", s, "library_leak", lib["loadings_exfiltrated"],
              "isolated_bypass_blocked", iso.get("natural_egress_blocked"),
              "native_blocked", native.get("native_egress_blocked"),
              "sanctioned_blocked", iso.get("sanctioned_undeclared_blocked"),
              "chain", transp.get("chain_verified"), flush=True)
    agg = dict(
        library_leak_frac=float(np.mean([r["library_arm"]["loadings_exfiltrated"] for r in per_seed])),
        natural_egress_blocked_frac=float(np.mean([r["isolated"]["natural_egress_blocked"] for r in per_seed])),
        native_egress_blocked_frac=float(np.mean([r["native"]["native_egress_blocked"] for r in per_seed])),
        sanctioned_blocked_frac=float(np.mean([r["isolated"]["sanctioned_undeclared_blocked"] for r in per_seed])),
        declared_flow_frac=float(np.mean([r["transport"]["declared_ok"] for r in per_seed])),
        chain_verified_frac=float(np.mean([r["transport"]["chain_verified"] for r in per_seed])),
        governed_af_auc=float(np.mean([r["transport"]["governed_af_auc"] for r in per_seed])),
        governed_pc_auc=float(np.mean([r["transport"]["governed_pc_auc"] for r in per_seed])),
        pergroup_full_delivered_frac=float(np.mean([r["pergroup"]["full"]["delivered"] for r in per_seed])),
        pergroup_2dp_delivered_frac=float(np.mean([r["pergroup"]["2dp"]["delivered"] for r in per_seed])),
        pergroup_chained_frac=float(np.mean([r["pergroup"]["full"]["chained"] and r["pergroup"]["2dp"]["chained"] for r in per_seed])),
        udp_egress_blocked_frac=float(np.mean([r["isolated"]["udp_egress_blocked"] is True for r in per_seed])))
    res = dict(per_seed=per_seed, agg=agg, n_seeds=len(seeds))
    with open(os.path.join(OUT, "s39_enforcement.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("agg", json.dumps({k: round(v, 3) for k, v in agg.items()}), flush=True)
    return res


if __name__ == "__main__":
    run()
