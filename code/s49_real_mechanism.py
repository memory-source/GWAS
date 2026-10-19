"""S49 - real within-population mechanism + granularity.

The paper's mechanism (the leak is the in-sample-ness of the study-fitted
spectrum; loadings not eigenvalues) is verified on real data only at the sub-AF
magnitude. s48 shows the real within-population PC channel is the dominant one
(gap up to +0.13). This re-runs the granularity + mechanism ladder on that
operating point:

  - top-k loadings (k=1,3,5,10) of the study-fitted basis;
  - eigenvalue-only control (random orthonormal basis, same k);
  - reference-fitted-basis control (an independent panel's loadings);
  - out-of-sample transfer (a held-out study half's basis applied to the other
    half, versus the same-half basis).

Cohorts: GBR (EUR) and CHB (EAS), chr22 + chr1. 5 split seeds.

Output: results/raw/s49_real_mechanism.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import allele_freq, af_lrt_attack, auroc, pca
from s48_real_structure import pop_table, load_for, _orient_auroc
from real_multi import SUPERPOP5

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")


def _pd(X, L, mean, M):
    return np.linalg.norm((X.astype(float) - mean) @ L.T / np.sqrt(M), axis=1)


def cell(G, lab, pop, seed, k=10):
    m = lab == pop
    idx = np.random.default_rng(seed).permutation(np.where(m)[0])
    q = idx.size // 3
    si, ri, ni = idx[:q], idx[q:2 * q], idx[2 * q:3 * q]
    site, ref, non = G[si], G[ri], G[ni]
    M = G.shape[1]
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    # one shared calibration/evaluation split for the main block; every channel
    # here is scored on the SAME held-out half, so PC-AF is a paired comparison.
    cal = np.random.default_rng(seed + 100).random(len(y)) < 0.5
    af_s, af_r = allele_freq(site), allele_freq(ref)
    am = np.concatenate([af_lrt_attack(site, af_s, af_r),
                         af_lrt_attack(non, af_s, af_r)])
    af_auc = _orient_auroc(am, y, np.random.default_rng(seed + 100), cal_mask=cal)
    _, _, L = pca(site, k=k); mean = site.mean(0)
    out = dict(af_auc=af_auc, n_study=int(len(site)), M=int(M))
    for K in (1, 3, 5, 10):
        dm, dn = _pd(site, L[:K], mean, M), _pd(non, L[:K], mean, M)
        out[f"top{K}"] = _orient_auroc(np.concatenate([dm, dn]), y,
                                       np.random.default_rng(seed + 100), cal_mask=cal)
    dm, dn = _pd(site, L, mean, M), _pd(non, L, mean, M)
    out["full"] = _orient_auroc(np.concatenate([dm, dn]), y,
                                np.random.default_rng(seed + 100), cal_mask=cal)
    out["gap"] = out["full"] - af_auc
    # eigenvalue-only
    rng = np.random.default_rng(seed + 5)
    Le = rng.normal(size=(k, M)); Le, _ = np.linalg.qr(Le.T); Le = Le.T
    dm, dn = _pd(site, Le, mean, M), _pd(non, Le, mean, M)
    out["eig_only"] = _orient_auroc(np.concatenate([dm, dn]), y,
                                    np.random.default_rng(seed + 100), cal_mask=cal)
    # reference-fitted basis
    _, _, Lr = pca(ref, k=k); mr = ref.mean(0)
    dm, dn = _pd(site, Lr, mr, M), _pd(non, Lr, mr, M)
    out["refbasis"] = _orient_auroc(np.concatenate([dm, dn]), y,
                                    np.random.default_rng(seed + 100), cal_mask=cal)
    # out-of-sample transfer: fit on one study half, apply to the OTHER half +
    # non-members (rows with no overlap between fit and score); a same-half basis
    # is the in-sample reference. Orientation uses a disjoint calibration half of
    # the scored rows.
    rng2 = np.random.default_rng(seed + 11)
    half = rng2.random(len(site)) < 0.5
    _, _, L_h = pca(site[half], k=k); m_h = site[half].mean(0)
    dm, dn = _pd(site[~half], L_h, m_h, M), _pd(non, L_h, m_h, M)
    y2 = np.concatenate([np.ones((~half).sum()), np.zeros(len(non))])
    cal2 = np.random.default_rng(seed + 140).random(len(y2)) < 0.5
    out["oos_transfer"] = _orient_auroc(np.concatenate([dm, dn]), y2,
                                        np.random.default_rng(seed + 140), cal_mask=cal2)
    _, _, L_s = pca(site[~half], k=k); m_s = site[~half].mean(0)
    dm, dn = _pd(site[~half], L_s, m_s, M), _pd(non, L_s, m_s, M)
    out["samehalf"] = _orient_auroc(np.concatenate([dm, dn]), y2,
                                    np.random.default_rng(seed + 140), cal_mask=cal2)
    return out


def run(pops=("GBR", "CHB"), chroms=(22,), n_splits=5):
    res = {}
    for chrom in chroms:
        G, lab, pos = load_for(SUPERPOP5, chrom, ld_window=1500, max_snps=20000)
        for pop in pops:
            rows = [cell(G, lab, pop, s) for s in range(n_splits)]
            agg = {}
            for k in rows[0]:
                if not isinstance(rows[0][k], (int, float)):
                    continue
                vals = np.array([r[k] for r in rows], float)
                agg[k] = float(vals.mean())
                if k not in ("n_study", "M"):
                    agg[k + "_std"] = float(vals.std())
            res[f"chr{chrom}:{pop}"] = agg
            print("%-12s AF %.3f full %.3f gap %+.3f top1/top3 %.3f/%.3f eig %.3f "
                  "refb %.3f oos %.3f same %0.3f" % (
                      f"chr{chrom}:{pop}", agg["af_auc"], agg["full"], agg["gap"],
                      agg["top1"], agg["top3"], agg["eig_only"], agg["refbasis"],
                      agg["oos_transfer"], agg["samehalf"]), flush=True)
    with open(os.path.join(OUT, "s49_real_mechanism.json"), "w") as f:
        json.dump(dict(cells=res, n_splits=n_splits, pops=list(pops),
                       chroms=list(chroms)), f, indent=2)
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default="22")
    ap.add_argument("--pops", default="GBR,CHB")
    ap.add_argument("--splits", type=int, default=5)
    a = ap.parse_args()
    run(tuple(a.pops.split(",")), tuple(int(c) for c in a.chroms.split(",")), a.splits)
