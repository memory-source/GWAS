"""S66 - control-provenance check for the no-release control (deep round).

The no-release control asks whether an attacker holding the known calibration members'
full genotypes (but no release) can separate members. In s61 the calibration PCA spans
are fitted and scored in the pipeline's public standardized space. This script tests
whether the control's conclusion depends on that preprocessing choice by running the
calibration-span control end to end in THREE coordinate systems. Each space transforms
the SPAN-FITTING data and the SCORED data with the SAME transform:

  public    x = (g - mu_pub) / sd_pub, with mu_pub / sd_pub the honest+partner public
            baseline (the pipeline's own preprocessing; matches s61 exactly)
  owncs     x = (g - mean_cal) / sd_cal, the calibration members' own per-marker mean
            and standard deviation (computed on the calibration members only)
  raw       x = g - mean_pool, no variance standardization (mean of the scored rows)

For each space: fit the top-k (and top-qk) PCA span of the calibration members'
transformed genotypes, score the pooled (members + non-members) rows' projection norm,
held-out oriented, and report the AUROC. |z| is bounded at s61.Z_CLIP in the two
standardized spaces; the raw space is unbounded.

Output: results/raw/s66_control_provenance_chr<C>.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
import s51_real_lowstruct as s51
import s61_spectral_tomography as s61
from common import auroc
from real_multi import SUPERPOP

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")


def cohorts_for(chrom, ld_window=1500, max_snps=20000):
    Ge, _, _ = s48.load_for(("EUR",), chrom, ld_window, max_snps)
    G3, _, sup3, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    G5, _, sup5, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    return {"within:EUR": (Ge, np.array(["EUR"] * len(Ge))),
            "super:3pop": (G3, sup3),
            "within:EAS": (G5[sup5 == "EAS"], sup5[sup5 == "EAS"])}


def _auc_from_span(R, X, y, cal, z_clip):
    """AUROC of the projection norm onto span(R) for rows X already in R's space,
    held-out oriented."""
    U, s, _ = np.linalg.svd(R, full_matrices=False)
    keep = s > 1e-9 * (s[0] if s.size else 1.0)
    B = U[:, keep]
    Z = np.clip(X, -z_clip, z_clip) if z_clip is not None else X
    d = np.linalg.norm(Z @ B, axis=1)
    h = np.asarray(cal, bool)
    sgn = 1.0 if auroc(d[h], y[h]) >= 0.5 else -1.0
    return float(auroc(sgn * d[~h], y[~h]))


def _run_space(Gcal, Gpool, y, cal, d_vals, space, k, seed, mu_pub, sd_pub,
               clip=s61.Z_CLIP):
    """Fit the calibration span(s) and score the pooled rows in ONE self-consistent
    space. Returns {d: auc} for each d in d_vals."""
    if space == "public":
        Tcal = (Gcal - mu_pub[None, :]) / sd_pub[None, :]
        Tpool = (Gpool - mu_pub[None, :]) / sd_pub[None, :]
        cl = clip
    elif space == "owncs":
        mu_c = Gcal.mean(0)
        sd_c = Gcal.std(0); sd_c[sd_c == 0] = 1.0
        Tcal = (Gcal - mu_c[None, :]) / sd_c[None, :]
        Tpool = (Gpool - mu_c[None, :]) / sd_c[None, :]
        cl = clip
    elif space == "raw":
        mu_p = Gpool.mean(0)
        Tcal = Gcal - mu_p[None, :]
        Tpool = Gpool - mu_p[None, :]
        cl = None
    else:
        raise ValueError(space)
    out = {}
    for d in d_vals:
        _, Vk = s61.topk([Tcal], d)
        out[d] = _auc_from_span(Vk, Tpool, y, cal, cl)
    return out


def one_split(G, groups, seed, k=s61.K_TOP, caps=None):
    roles = s61.split_roles(groups, caps or s61.CAPS, seed)
    honest, partner = G[roles["honest"]], G[roles["partner"]]
    nonmem = G[roles["nonmember"]]
    mu, sd = s61.joint_baseline(honest, partner)          # the public baseline
    y = np.concatenate([np.ones(len(honest)), np.zeros(len(nonmem))])
    cal = np.random.default_rng(seed + 31).random(len(y)) < 0.5
    h = np.asarray(cal, bool)
    # calibration members = the honest rows in the calibration half ONLY (matches the
    # s61 control: `Zk = Zh[h[:len(honest)]]`, i.e. honest calibration rows)
    Gcal = honest[h[:len(honest)]].astype(np.float64)
    Gpool = np.vstack([honest, nonmem]).astype(np.float64)
    d_k = k
    d_qk = min(max(s61.QS) * k, Gcal.shape[0] - 1)
    out = {}
    for space in ("public", "owncs", "raw"):
        res = _run_space(Gcal, Gpool, y, cal, (d_k, d_qk), space, k, seed, mu, sd)
        out[f"cal_pca_k_{space}"] = res[d_k]
        out[f"cal_pca_qk_{space}"] = res[d_qk]
        # random-span control in the same space (row transform only sets the clip)
        rng = np.random.default_rng(seed + 71)
        Rrand = rng.standard_normal((Gcal.shape[1], k))
        if space == "public":
            Tpool = (Gpool - mu[None, :]) / sd[None, :]; cl = s61.Z_CLIP
        elif space == "owncs":
            mu_c = Gcal.mean(0); sd_c = Gcal.std(0); sd_c[sd_c == 0] = 1.0
            Tpool = (Gpool - mu_c[None, :]) / sd_c[None, :]; cl = s61.Z_CLIP
        else:
            Tpool = Gpool - Gpool.mean(0)[None, :]; cl = None
        out[f"random_k_{space}"] = _auc_from_span(Rrand, Tpool, y, cal, cl)
    return out


def run(chrom=22, seeds=5):
    out = {"chrom": chrom, "seeds": seeds, "cohorts": {}}
    for name, (G, groups) in cohorts_for(chrom).items():
        rows = [one_split(G, groups, s) for s in range(seeds)]
        keys = rows[0].keys()
        blk = {k: float(np.mean([r[k] for r in rows])) for k in keys}
        for k in keys:
            blk[k + "_ci95"] = s61._boot_ci(np.array([r[k] for r in rows]))
        out["cohorts"][name] = blk
        print("### %s" % name, flush=True)
        for k in keys:
            print("   %-24s %.3f" % (k, blk[k]), flush=True)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s66_control_provenance_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
