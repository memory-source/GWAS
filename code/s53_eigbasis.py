"""S53 - fixed (non-random) eigenvalue-informed projection control.

The paper's eigenvalue-only negative is currently a RANDOM-basis projection: an
attacker who knows only the top-k eigenvalues substitutes a random orthonormal
basis. A reviewer could object that the negative is an artifact of randomness.
This probes a fixed, DATA-INDEPENDENT basis that DOES carry the study's
eigenvalue spectrum:

  - `lam_weight` basis: the first-k standard axes (a fixed, exactly known
    orthonormal basis) weighted by the released eigenvalue spectrum, so the
    projection distance `d(x) = sum_j ((x-mu).b_j)^2 / lambda_j` is a
    fixed-basis, eigenvalue-inspired Mahalanobis-like score. This uses the
    released eigenvalues (which the random-basis control discards) but no
    study-fitted loading directions.

Plus two zero-information fixed bases for reference:

  - the first-k standard-axis (identity) basis, unweighted;
  - a deterministic cosine-direction basis (a fixed, smooth, data-independent set
    of directions).

All are scored like the eigenvalue-only control (projection distance to the study
centroid). If the eigenvalue-weighted fixed basis is at chance, the eigenvalue-only
negative does not depend on the basis being random AND is not rescued by reusing
the released eigenvalue magnitudes. Full study-fitted loadings are the positive
reference.

Cells: a simulated 3-group panel and pooled real within-super-population cells
(chr22/chr1/chr2 x EUR/EAS/AFR). 5 seeds.

Output: results/raw/s53_eigbasis.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc, pca)
from s48_real_structure import load_for, _orient_auroc
from real_multi import SUPERPOP5

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")


def _pd(X, L, mean):
    Xf = X.astype(np.float64) - mean[None, :]
    return np.linalg.norm((Xf @ L.T) / np.sqrt(max(1, X.shape[1])), axis=1)


def _lam_weight_pd(X, mean, k, lam):
    """Fixed-basis, eigenvalue-weighted projection distance:
    d(x) = sum_{j<k} ((x-mu).b_j)^2 / lambda_j, b_j the standard axis j.
    Uses the released eigenvalues, no study-fitted directions."""
    Xf = X.astype(np.float64) - mean[None, :]
    lam = np.clip(np.asarray(lam, float)[:k], 1e-9, None)
    # standard axes: component j is Xf[:, j]; weighted Mahalanobis-like distance
    W = np.zeros(Xf.shape[1])
    W[:k] = 1.0 / lam
    return np.sqrt((Xf ** 2 @ W))


def _dct_basis(M, k):
    """First-k deterministic cosine-direction basis vectors, (k, M): fixed and
    data-independent (a smooth cosine set, not exactly an orthonormal DCT)."""
    n = np.arange(M)
    B = np.empty((k, M))
    for i in range(k):
        B[i] = np.cos(np.pi * (i + 0.5) * n / M)
    B /= np.linalg.norm(B, axis=1, keepdims=True)
    return B


def _identity_basis(M, k):
    B = np.zeros((k, M))
    for i in range(k):
        B[i, i] = 1.0
    return B


def cell_sim(seed, k=10):
    pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
    site, ref, non = split_population3(pop, 300, 300, 300)
    rel = release_objects(site)
    mask = rel["mask"]
    G = site["G"][:, mask]
    Gn = non["G"][:, mask]
    M = mask.sum()
    y = np.concatenate([np.ones(len(G)), np.zeros(len(Gn))])
    cal = np.random.default_rng(seed + 100).random(len(y)) < 0.5
    eig, _, L = pca(site["G"], k=k)
    L = L[:, mask]
    mean = G.mean(0)
    out = dict(seed=seed, M=int(M))
    out["full"] = _orient_auroc(np.concatenate([_pd(G, L, mean), _pd(Gn, L, mean)]),
                                y, np.random.default_rng(seed + 100), cal_mask=cal)
    rng = np.random.default_rng(seed + 5)
    Le = rng.normal(size=(k, M)); Le, _ = np.linalg.qr(Le.T); Le = Le.T
    out["rand_basis"] = _orient_auroc(np.concatenate([_pd(G, Le, mean), _pd(Gn, Le, mean)]),
                                      y, np.random.default_rng(seed + 100), cal_mask=cal)
    Bi = _identity_basis(M, k)
    out["identity_basis"] = _orient_auroc(np.concatenate([_pd(G, Bi, mean), _pd(Gn, Bi, mean)]),
                                          y, np.random.default_rng(seed + 100), cal_mask=cal)
    Bd = _dct_basis(M, k)
    out["dct_basis"] = _orient_auroc(np.concatenate([_pd(G, Bd, mean), _pd(Gn, Bd, mean)]),
                                     y, np.random.default_rng(seed + 100), cal_mask=cal)
    lw = np.concatenate([_lam_weight_pd(G, mean, k, eig), _lam_weight_pd(Gn, mean, k, eig)])
    out["lam_weight_basis"] = _orient_auroc(lw, y, np.random.default_rng(seed + 100), cal_mask=cal)
    return out


def cell_real(G, crit, sp, seed, k=10):
    m = crit == sp
    Gp = G[m]
    idx = np.random.default_rng(seed).permutation(len(Gp))
    q = idx.size // 3
    if q < 6:
        return None
    si, ri, ni = idx[:q], idx[q:2 * q], idx[2 * q:3 * q]
    site, ref, non = Gp[si], Gp[ri], Gp[ni]
    M = G.shape[1]
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    cal = np.random.default_rng(seed + 100).random(len(y)) < 0.5
    eig, _, L = pca(site, k=k)
    mean = site.mean(0)
    out = dict(seed=seed, M=int(M), n_study=int(len(site)))
    out["full"] = _orient_auroc(np.concatenate([_pd(site, L, mean), _pd(non, L, mean)]),
                                y, np.random.default_rng(seed + 100), cal_mask=cal)
    rng = np.random.default_rng(seed + 5)
    Le = rng.normal(size=(k, M)); Le, _ = np.linalg.qr(Le.T); Le = Le.T
    out["rand_basis"] = _orient_auroc(np.concatenate([_pd(site, Le, mean), _pd(non, Le, mean)]),
                                      y, np.random.default_rng(seed + 100), cal_mask=cal)
    Bi = _identity_basis(M, k)
    out["identity_basis"] = _orient_auroc(np.concatenate([_pd(site, Bi, mean), _pd(non, Bi, mean)]),
                                          y, np.random.default_rng(seed + 100), cal_mask=cal)
    Bd = _dct_basis(M, k)
    out["dct_basis"] = _orient_auroc(np.concatenate([_pd(site, Bd, mean), _pd(non, Bd, mean)]),
                                     y, np.random.default_rng(seed + 100), cal_mask=cal)
    lw = np.concatenate([_lam_weight_pd(site, mean, k, eig), _lam_weight_pd(non, mean, k, eig)])
    out["lam_weight_basis"] = _orient_auroc(lw, y, np.random.default_rng(seed + 100), cal_mask=cal)
    return out


def _agg(rows, skip=("seed", "M", "n_study")):
    out = {}
    for k in rows[0]:
        if k in skip or not isinstance(rows[0][k], (int, float)):
            continue
        vals = np.array([r[k] for r in rows], float)
        out[k] = float(vals.mean())
        out[k + "_std"] = float(vals.std())
    return out


def run(n_splits=5, chroms=(22, 1, 2), sups=("EUR", "EAS", "AFR")):
    res = {}
    rows = [cell_sim(s) for s in range(n_splits)]
    res["sim_3group"] = _agg(rows)
    print("sim_3group  full %.3f  rand %.3f  identity %.3f  dct %.3f  lamW %.3f" % (
        res["sim_3group"]["full"], res["sim_3group"]["rand_basis"],
        res["sim_3group"]["identity_basis"], res["sim_3group"]["dct_basis"],
        res["sim_3group"]["lam_weight_basis"]), flush=True)
    samples, sups_all, pops = __import__("s48_real_structure").pop_table()
    order = [s for s, kk in zip(samples, sups_all) if kk in set(SUPERPOP5)]
    crit = np.array([sups_all[samples.index(s)] for s in order])
    for chrom in chroms:
        G, lab, pos = load_for(SUPERPOP5, chrom, ld_window=1500, max_snps=20000)
        for sp in sups:
            rr = [cell_real(G, crit, sp, s) for s in range(n_splits)]
            rr = [r for r in rr if r]
            if rr:
                res[f"real_chr{chrom}:{sp}"] = _agg(rr)
                a = res[f"real_chr{chrom}:{sp}"]
                print("real_chr%d:%-4s full %.3f  rand %.3f  identity %.3f  dct %.3f  lamW %.3f" % (
                    chrom, sp, a["full"], a["rand_basis"], a["identity_basis"], a["dct_basis"],
                    a["lam_weight_basis"]), flush=True)
    with open(os.path.join(OUT, "s53_eigbasis.json"), "w") as f:
        json.dump(dict(cells=res, n_splits=n_splits, chroms=list(chroms),
                       sups=list(sups)), f, indent=2)
    return res


if __name__ == "__main__":
    run()
