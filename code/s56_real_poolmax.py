"""S56 - real structure-zero switch, balanced design (round-16 B1).

The pooled single-super-population arm of s51 holds measured between-group
structure exactly zero (one group) but caps the study at n_study=167-220 under an
equal thirds split. This grows the study toward the pool maximum while HOLDING
THE REFERENCE PANEL SIZE FIXED, so the AF-LRT baseline is not disadvantaged by a
shrinking reference (a study-heavy split that also shrank the reference would
lower AF and inflate the PC-AF gap mechanically). Non-members absorb the
remainder.

Schedule (n_study, n_ref); non = pool - n_study - n_ref:
  (167, 150), (250, 150), (300, 150), (350, 150), and the pool-max study with the
  reference still 150 where the pool allows (AFR). Controls: reference-fitted
  basis and eigenvalue-only. 5 split seeds; per-split gaps exported.

Output: results/raw/s56_real_poolmax_chr<C>.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
from common import allele_freq, af_lrt_attack, auroc, pca
from real_multi import SUPERPOP5

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
N_SPLITS = 5
N_REF = 150
STUDY_SCHEDULE = (167, 250, 300, 350)


def _orient(scores, y, cal_mask):
    return s48._orient_auroc(scores, y, np.random.default_rng(0), cal_mask=cal_mask)


def cell(G, n_study, seed, k=10):
    rng = np.random.default_rng(seed)
    n = G.shape[0]
    # reference held FIXED at N_REF, non-members absorb the remainder; where the
    # pool cannot supply N_REF non-members the reference is capped (its size is
    # recorded per cell) so the AF baseline is not shrunk by design.
    nr = min(N_REF, n - n_study - 6)
    if nr < 6 or n_study < 30:
        return None
    idx = rng.permutation(n)
    si = idx[:n_study]; ri = idx[n_study:n_study + nr]; ni = idx[n_study + nr:]
    if len(ni) < 6:
        return None
    site, ref, non = G[si], G[ri], G[ni]
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    cal = np.random.default_rng(seed + 31).random(len(y)) < 0.5
    af_s, af_r = allele_freq(site), allele_freq(ref)
    am = np.concatenate([af_lrt_attack(site, af_s, af_r), af_lrt_attack(non, af_s, af_r)])
    af_auc = _orient(am, y, cal)
    _, _, L = pca(site, k=k); mean = site.mean(0)
    dm = np.linalg.norm((site.astype(float) - mean) @ L.T / np.sqrt(G.shape[1]), axis=1)
    dn = np.linalg.norm((non.astype(float) - mean) @ L.T / np.sqrt(G.shape[1]), axis=1)
    pc_auc = _orient(np.concatenate([dm, dn]), y, cal)
    _, _, Lr = pca(ref, k=k); mr = ref.mean(0)
    xm = np.linalg.norm((site.astype(float) - mr) @ Lr.T / np.sqrt(G.shape[1]), axis=1)
    xn = np.linalg.norm((non.astype(float) - mr) @ Lr.T / np.sqrt(G.shape[1]), axis=1)
    refb = _orient(np.concatenate([xm, xn]), y, cal)
    rng2 = np.random.default_rng(seed + 5)
    Le = rng2.normal(size=(k, G.shape[1])); Le, _ = np.linalg.qr(Le.T); Le = Le.T
    em = np.linalg.norm((site.astype(float) - mean) @ Le.T / np.sqrt(G.shape[1]), axis=1)
    en = np.linalg.norm((non.astype(float) - mean) @ Le.T / np.sqrt(G.shape[1]), axis=1)
    eigo = _orient(np.concatenate([em, en]), y, cal)
    return dict(af_auc=af_auc, pc_auc=pc_auc, gap=pc_auc - af_auc,
                pc_refbasis=refb, pc_eigonly=eigo,
                n_study=int(len(site)), n_ref=int(len(ref)), n_non=int(len(non)),
                M=int(G.shape[1]))


def _agg(rows):
    out = {}
    for k in rows[0]:
        if not isinstance(rows[0][k], (int, float)):
            continue
        vals = np.array([r[k] for r in rows], float)
        out[k] = float(vals.mean())
        if k not in ("n_study", "n_ref", "n_non", "M"):
            out[k + "_std"] = float(vals.std())
    return out


def run(chrom=22, sps=("EUR", "EAS", "AFR")):
    res = {"chrom": chrom, "n_splits": N_SPLITS, "n_ref_target": N_REF, "cells": {}}
    G5, lab5, _ = s48.load_for(SUPERPOP5, chrom, ld_window=1500, max_snps=20000)
    samples, sups, pops = s48.pop_table()
    order = [s for s, kk in zip(samples, sups) if kk in set(SUPERPOP5)]
    order_sup = np.array([sups[samples.index(s)] for s in order])
    for sp in sps:
        m = order_sup == sp
        n_pool = int(m.sum())
        for ns in STUDY_SCHEDULE:
            if ns + 30 > n_pool:
                continue
            rows = [cell(G5[m], ns, s) for s in range(N_SPLITS)]
            rows = [r for r in rows if r]
            if not rows:
                continue
            key = f"{sp}:n{ns}"
            a = _agg(rows)
            a["n_pool"] = n_pool
            a["per_split_gap"] = [r["gap"] for r in rows]
            res["cells"][key] = a
            print("%-12s n_pool %3d n_study %3d n_ref %3d  AF %.3f  PC %.3f  "
                  "gap %+.3f (sd %.3f)  refb %.3f  eigo %.3f" % (
                      key, n_pool, a["n_study"], a["n_ref"], a["af_auc"], a["pc_auc"],
                      a["gap"], a["gap_std"], a["pc_refbasis"], a["pc_eigonly"]), flush=True)
            with open(os.path.join(OUT, f"s56_real_poolmax_chr{chrom}.json"), "w") as f:
                json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s56_real_poolmax_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
