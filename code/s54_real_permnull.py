"""S54 - permutation nulls and per-split dispersion for the REAL-data gaps.

Round-16/17 flagged the one cheap strengthening left: the real round-13/15 gaps
(s48 within-pop, s48 3/5-pop, s50 3/5-pop, s51 pooled and fixed-composition
ladder) are reported as five-split means with no dispersion and no null. This
adds, on the SAME cohorts and the SAME held-out orientation as s48
(`cohort_cell` / `_orient_auroc`), for every cell:

  - the five per-split gaps, their std, and a bootstrap CI with P(gap > 0);
  - a stratified label-permutation null on the *oriented evaluation half* (the
    sign is fixed on the disjoint calibration half exactly as in s48, and labels
    are permuted within the cohort's group among the evaluated rows), giving an
    empirical p-value that each channel separates members.

AF and PC use the identical attack code path as s48, so the observed AUROC and
the gap reproduce the s48 JSONs; the s48 JSONs themselves are untouched.

Output: results/raw/s54_real_permnull.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
from common import allele_freq, af_lrt_attack, auroc, pca
from real_multi import SUPERPOP, SUPERPOP5

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
N_SPLITS = 5
N_PERM = 5000


def _oriented(scores, y, cal_mask):
    """Return (oriented_eval_scores, y_eval, strata_ok). Mirrors s48._orient_auroc:
    the sign is chosen on the calibration half and the AUROC evaluated on the
    other, disjoint half; no scored individual's label selects the direction."""
    s = np.asarray(scores, float)
    h = np.asarray(cal_mask, bool)
    if h.sum() >= 2 and (~h).sum() >= 2 and len(np.unique(y[h])) == 2:
        sgn = 1.0 if auroc(s[h], y[h]) >= 0.5 else -1.0
        ev = ~h
    else:
        sgn = -1.0
        ev = np.ones(len(s), bool)
    return sgn * s[ev], np.asarray(y)[ev]


def _perm_pvalue(ev_scores, y_ev, strata_ev, n_perm, seed):
    """P(permuted AUROC >= observed), permuting labels within strata on the
    evaluated rows only (scores fixed)."""
    obs = auroc(ev_scores, y_ev)
    rng = np.random.default_rng(seed)
    groups = [np.where(strata_ev == g)[0] for g in np.unique(strata_ev)]
    ge = 0
    for _ in range(n_perm):
        perm = y_ev.copy()
        for idx in groups:
            perm[idx] = rng.permutation(perm[idx])
        ge += int(auroc(ev_scores, perm) >= obs)
    return float(obs), (ge + 1) / (n_perm + 1)


def _cell(G, lab, per, seeds, k=10):
    afs, pcs, af_p, pc_p = [], [], [], []
    for s in range(seeds):
        rng = np.random.default_rng(s)
        si, ri, ni = [], [], []
        for g in np.unique(lab):
            idx = rng.permutation(np.where(lab == g)[0])
            q = min(per, idx.size // 3)
            if q < 4:
                continue
            si.append(idx[:q]); ri.append(idx[q:2 * q]); ni.append(idx[2 * q:3 * q])
        if not si:
            return None
        si, ri, ni = np.concatenate(si), np.concatenate(ri), np.concatenate(ni)
        site, ref, non = G[si], G[ri], G[ni]
        y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
        strata = np.concatenate([lab[si], lab[ni]])
        cal = np.random.default_rng(s + 31).random(len(y)) < 0.5
        af_s, af_r = allele_freq(site), allele_freq(ref)
        af_raw = np.concatenate([af_lrt_attack(site, af_s, af_r),
                                 af_lrt_attack(non, af_s, af_r)])
        _, _, L = pca(site, k=k); mean = site.mean(0)
        pc_raw = np.linalg.norm(
            (np.vstack([site, non]).astype(float) - mean) @ L.T / np.sqrt(G.shape[1]), axis=1)
        afe, aye = _oriented(af_raw, y, cal)
        pce, _ = _oriented(pc_raw, y, cal)
        afo, afp = _perm_pvalue(afe, aye, strata[~cal], N_PERM, seed=s + 1)
        pco, pcp = _perm_pvalue(pce, aye, strata[~cal], N_PERM, seed=s + 2)
        afs.append(afo); pcs.append(pco); af_p.append(afp); pc_p.append(pcp)
    afs, pcs = np.array(afs), np.array(pcs)
    gaps = pcs - afs
    rng = np.random.default_rng(123)
    boot = gaps[rng.integers(0, len(gaps), size=(10000, len(gaps)))].mean(1)
    return dict(af_auc=float(afs.mean()), pc_auc=float(pcs.mean()),
                gap=float(gaps.mean()), gap_std=float(gaps.std()),
                gap_ci=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                gap_p_gt0=float((gaps > 0).mean()),
                af_p_max=float(max(af_p)), pc_p_max=float(max(pc_p)),
                n_splits=int(len(gaps)))


def run(chrom=22, n_splits=N_SPLITS):
    res = {"chrom": chrom, "n_splits": n_splits, "n_perm": N_PERM, "cells": {}}
    G5, lab5, _ = s48.load_for(SUPERPOP5, chrom, ld_window=1500, max_snps=20000)
    samples, sups, pops = s48.pop_table()
    order = [s for s, kk in zip(samples, sups) if kk in set(SUPERPOP5)]
    order_sup = np.array([sups[samples.index(s)] for s in order])
    # per=150 matches s48.cohort_cell's default so the observed gaps reproduce
    # the s48 JSONs; for a single population the cap (size//3) binds below 150.
    for pop in ("GBR", "CHB", "YRI"):
        m = lab5 == pop
        if m.sum() >= 12:
            r = _cell(G5[m], lab5[m], per=150, seeds=n_splits)
            if r:
                res["cells"][f"within:{pop}"] = r
    for name, sp_keep in (("3pop", SUPERPOP), ("5pop", SUPERPOP5)):
        m = np.isin(order_sup, sp_keep)
        r = _cell(G5[m], order_sup[m], per=150, seeds=n_splits)
        if r:
            res["cells"][name] = r
    with open(os.path.join(OUT, "s54_real_permnull.json"), "w") as f:
        json.dump(res, f, indent=2)
    for name, a in res["cells"].items():
        print("%-14s AF %.3f(p~%.0e) PC %.3f(p~%.0e)  gap %+.3f (sd %.3f) CI[%+.3f,%+.3f] P>0 %.2f" % (
            name, a["af_auc"], a["af_p_max"], a["pc_auc"], a["pc_p_max"],
            a["gap"], a["gap_std"], a["gap_ci"][0], a["gap_ci"][1], a["gap_p_gt0"]),
            flush=True)
    return res


if __name__ == "__main__":
    t0 = time.time()
    payload = run()
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "s54_real_permnull.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
