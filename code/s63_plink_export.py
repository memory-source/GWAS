"""S63 - real-software export interface check via PLINK 2.0 --pca.

Runs the selection-input spectral tomography through a REAL PCA export. PLINK 2.0
extracts the joint PCs of (honest, partner_t) and writes the allele weights
(`--pca <k> allele-wts` -> `.eigenvec.allele`, one row per allele) and eigenvalues
(`.eigenval`). The attacker then READS BACK those files and forms, in PLINK's own
variance-standardized scaling,

    C_H = Z_H^T Z_H / M        (Z standardized as in PLINK: (g-2p)/sqrt(2p(1-p)))
    K_t = Z_P(t)^T Z_P(t) / M
    Y_t = V_t Lambda_t - K_t V_t   =   C_H V_t

PLINK's eigenvalues equal the scatter eigenvalues divided by M, so the identity
holds in PLINK's normalization. The AF-preserving block swap leaves p unchanged, so
the scaling is fixed across releases.

The script reports (i) the identity residual against the read-back files, (ii) the
marker order check, (iii) the membership AUROC from the stacked read-back outputs
versus the PLINK-conditioned single release, and (iv) the effect of the printed
6-significant-digit export precision.

Output: results/raw/s63_plink_export_chr<C>.json
"""
import json
import os
import subprocess
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
import s61_spectral_tomography as s61
from real_multi import load_panel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "raw")
PLINK = os.path.join(ROOT, "scratch", "plink2bin", "plink2")
WORK = os.path.join(ROOT, "scratch", "plink_export")
K = 10
N_BOOT = 10000


def write_vcf(path, G, pos, sample_ids):
    """Minimal VCF; dosages 0/1/2 -> GT 0/0,0/1,1/1. The ID column carries the
    marker position so the exported order can be checked."""
    gt = np.array([["0/0", "0/1", "1/1"][int(v)] for v in G.reshape(-1)]).reshape(G.shape)
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(sample_ids) + "\n")
        for j in range(G.shape[1]):
            f.write("22\t%d\t%s\tA\tG\t.\tPASS\t.\tGT\t" % (int(pos[j]), int(pos[j]))
                    + "\t".join(gt[:, j]) + "\n")


def run_plink_pca(vcf, out_prefix, k=K):
    cmd = [PLINK, "--vcf", vcf, "--allow-no-sex", "--pca", str(k), "allele-wts",
           "--out", out_prefix]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"plink2 failed rc={p.returncode}\n{p.stdout[-800:]}\n{p.stderr[-800:]}")


def read_pca(out_prefix):
    """Read the PLINK allele weights and eigenvalues. `.eigenvec.allele` columns are
    CHROM ID REF ALT A1 PC1..PCk with one row per allele (REF and ALT); the .pvar
    order equals the VCF order, and the ALT rows carry the alleles in panel order."""
    with open(out_prefix + ".eigenval") as f:
        lam = np.array([float(x) for x in f.read().split()])
    ids, rows = [], []
    with open(out_prefix + ".eigenvec.allele") as f:
        f.readline()
        for line in f:
            p = line.split()
            ids.append(p[1])
            rows.append([float(x) for x in p[5:]])
    rows = np.array(rows)
    V = rows[1::2]                                     # ALT rows, in marker order
    ids = ids[1::2]
    return lam[:K], V, ids


def std_panel(G, p, scale):
    return (G.astype(np.float64) - 2.0 * p[None, :]) / scale[None, :]


def membership_auc(R, honest, nonmem, mu2, scale, y, cal, sign=None):
    """AUROC and TPR@5%FPR from the range of a recovered operator R (marker side),
    held-out oriented. `sign` reuses an orientation fixed elsewhere. Returns
    (auc, tpr05, sign)."""
    U, s, _ = np.linalg.svd(R, full_matrices=False)
    keep = s > 1e-9 * (s[0] if s.size else 1.0)
    B = U[:, keep]
    X = np.vstack([honest, nonmem]).astype(np.float64)
    X = (X - 2.0 * mu2[None, :]) / scale[None, :]
    d = np.linalg.norm(X @ B, axis=1)
    h = np.asarray(cal, bool)
    sgn = sign if sign is not None else (1.0 if s61.auroc(d[h], y[h]) >= 0.5 else -1.0)
    s_ev, y_ev = sgn * d[~h], y[~h]
    return float(s61.auroc(s_ev, y_ev)), float(s61.tpr_at_fpr(s_ev, y_ev, 0.05)), sgn


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def run(chrom=22, seeds=3, block_len=None, caps=None):
    block_len = block_len or s61.BLOCK_LEN
    caps = caps or s61.CAPS
    G, lab, pos = s48.load_for(("EUR",), chrom, 1500, 20000)
    groups = np.array(["EUR"] * len(G))
    samples, sups = load_panel()
    sample_ids = [s for s, k in zip(samples, sups) if k == "EUR"]
    os.makedirs(WORK, exist_ok=True)
    n_blocks = G.shape[1] // block_len
    M = G.shape[1]

    identity, per_seed = [], []
    for s in range(seeds):
        roles = s61.split_roles(groups, caps, s)
        honest = G[roles["honest"]]; partner = G[roles["partner"]]
        nonmem = G[roles["nonmember"]]
        joint_ids = [sample_ids[i] for i in roles["honest"]] + \
                    [sample_ids[i] for i in roles["partner"]]
        p = np.vstack([honest, partner]).mean(0) / 2.0
        scale = np.sqrt(np.clip(2.0 * p * (1 - p), 1e-9, None))
        Zh = std_panel(honest, p, scale)
        C_H = Zh.T @ Zh / M
        y = np.concatenate([np.ones(len(honest)), np.zeros(len(nonmem))])
        cal = np.random.default_rng(s + 31).random(len(y)) < 0.5

        rng = np.random.default_rng(s + 101)
        Qs, Ys, lam0_V0 = [], [], None
        for t in range(8):                                 # outputs q = 1,2,4,8
            own = partner if t == 0 else s61.swap_block(
                partner, rng.choice(partner.shape[0], 2, replace=False),
                int(rng.integers(0, n_blocks)), block_len)[0]
            jp = np.vstack([honest, own])
            vcf = os.path.join(WORK, f"joint_s{s}_t{t}.vcf")
            write_vcf(vcf, jp, pos, joint_ids)
            pref = os.path.join(WORK, f"joint_s{s}_t{t}")
            run_plink_pca(vcf, pref)
            lam, V, ids_out = read_pca(pref)
            ok_order = all(ids_out[j] == str(int(pos[j])) for j in range(len(ids_out)))
            Zp = std_panel(own, p, scale)
            Kp = Zp.T @ Zp / M
            Y = V * lam[None, :] - Kp @ V
            residual = float(np.max(np.abs(C_H @ V - Y)) / max(1e-12, np.max(np.abs(C_H @ V))))
            identity.append({"seed": s, "q": t + 1, "marker_order_ok": bool(ok_order),
                             "rel_identity_residual": residual})
            Qs.append(V); Ys.append(Y)
            if t == 0:
                lam0_V0 = (lam, V, Kp)
        # membership on the stacked read-back outputs, and the strongest single
        # release (PLINK-conditioned: subtract the calibration members' scatter)
        aucs = {}
        tprs = {}
        for q in (1, 2, 4, 8):
            Q = np.concatenate(Qs[:q], axis=1); Y = np.concatenate(Ys[:q], axis=1)
            aucs[f"recovered_q{q}"], tprs[f"recovered_q{q}"], _ = membership_auc(
                Y @ np.linalg.pinv(Q), honest, nonmem, p, scale, y, cal)
        lam, V0, K0 = lam0_V0
        h = np.asarray(cal, bool)
        Zk = Zh[h[:len(honest)]]
        Y0 = V0 * lam[None, :] - K0 @ V0
        base, base_tpr, sgn0 = membership_auc(Y0, honest, nonmem, p, scale, y, cal)
        # the conditioning removes the calibration members, so reuse the sign fixed
        # on the unconditioned single release (as in s61/s10_paired)
        Y0c = Y0 - Zk.T @ (Zk @ V0) / M
        spec_cond, spec_cond_tpr, _ = membership_auc(Y0c, honest, nonmem, p, scale, y,
                                                     cal, sign=sgn0)
        # AF baselines under the same allowed information as s61/s62 (reference /
        # pooled background, and the composition-balanced conditioned AF)
        honest_groups = np.asarray(groups)[roles["honest"]]
        af_ref, af_pooled, af_cond, _nbal = s61.af_baselines(
            honest, nonmem, G[roles["reference"]], G[roles["partner"]],
            G[roles["shadow"]], honest_groups, y, cal, s)
        # strongest passive baseline: conditioned AF vs the conditioned single
        # release (both use the calibration members; the single release is the
        # spectral reader's own best passive channel)
        strong = max(spec_cond, af_cond, af_ref)
        per_seed.append({"seed": s, "af_auc": af_ref, "af_pooled": af_pooled,
                         "af_cond": af_cond, "single_release": base,
                         "single_release_tpr05": base_tpr,
                         "spec_cond_single": spec_cond,
                         "spec_cond_tpr05": spec_cond_tpr,
                         "tpr_fpr05": tprs, **aucs,
                         "gap_q8_vs_strongest_single": aucs["recovered_q8"] - strong})
        print("seed %d  af %.3f  afcond %.3f  single %.3f  spec_cond %.3f  q1 %.3f  q2 %.3f "
              "q4 %.3f  q8 %.3f  gap_q8 %+.3f" % (
                  s, af_ref, af_cond, base, spec_cond, aucs["recovered_q1"],
                  aucs["recovered_q2"], aucs["recovered_q4"], aucs["recovered_q8"],
                  aucs["recovered_q8"] - strong), flush=True)

    res = {"chrom": chrom, "seeds": seeds, "k": K,
           "plink_version": subprocess.run([PLINK, "--version"], capture_output=True,
                                           text=True).stdout.strip(),
           "normalization": "variance-standardized (g-2p)/sqrt(2p(1-p)); eigenvalues "
                            "are scatter-eigenvalues divided by marker count",
           "identity": identity, "per_seed": per_seed}
    gaps = np.array([r["gap_q8_vs_strongest_single"] for r in per_seed])
    rng = np.random.default_rng(7)
    boot = gaps[rng.integers(0, len(gaps), size=(N_BOOT, len(gaps)))].mean(1)
    res["gap_q8_mean"] = float(gaps.mean())
    res["gap_q8_ci95"] = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]
    res["gap_q8_pos_splits"] = int((gaps > 0).sum())
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    ap.add_argument("--seeds", type=int, default=3)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom, a.seeds)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s63_plink_export_chr{a.chrom}.json"), "w") as f:
        json.dump(_jsonable(payload), f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
