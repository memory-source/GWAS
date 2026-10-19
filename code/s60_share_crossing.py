"""S60 - cumulative squared-projection share for the s51 ladder + per-seed
crossing location.

The s51 ladder was quoted as a PC projection-distance gap; a |projection| share
can rise with k by construction. This recomputes, on the s51 fixed-composition
ladder (per-sample SUPER-POP label), the CUMULATIVE squared-projection share

      cum_k = sum_{j<=k} ||proj_j||^2 / sum_j ||proj_j||^2

(the running sum of the per-loading squared-projection contributions over the
total squared norm about the study mean), which is non-decreasing by construction
of a cumulative sum of non-negative increments. It also propagates the per-seed
crossing location of each ladder (linear interpolation between the two adjacent
rungs whose gap straddles zero) and reports its spread. Seeds whose gap sequence
crosses more than once (or never) are reported explicitly rather than silently
dropped.

Output: results/raw/s60_share_crossing.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
import s51_real_lowstruct as s51
from common import pca
from real_multi import SUPERPOP

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
N_SPLITS = 5
LADDER = (40, 80, 120, 160)
KM = 10


def _proj_norms(site, non, k):
    """Return squared-projection l2^2 of member and non-member rows onto the
    study-fitted top-k loadings (centered, unnormalized l2^2)."""
    _, _, L = pca(site, k=k)
    mean = site.mean(0)
    Pm = (site.astype(np.float64) - mean) @ L.T
    Pn = (non.astype(np.float64) - mean) @ L.T
    return Pm, Pn


def share_ladder(G, lab, per, seed):
    """Monotone per-k share of the squared projection for ONE split."""
    rng = np.random.default_rng(seed)
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
    site, non = G[si], G[ni]
    Pm, Pn = _proj_norms(site, non, KM)
    tot_m = (site.astype(np.float64) - site.mean(0)) ** 2
    tot_n = (non.astype(np.float64) - site.mean(0)) ** 2
    den_m = tot_m.sum(1); den_n = tot_n.sum(1)
    out = {}
    for k in range(1, KM + 1):
        inc_m = (Pm[:, k - 1] ** 2)
        inc_n = (Pn[:, k - 1] ** 2)
        # share of the k-th loading contribution, aggregated over individuals
        out[f"k{k}"] = float((inc_m.sum() + inc_n.sum()) /
                             (den_m.sum() + den_n.sum() + 1e-12))
    # true monotone share; the cumulative is what makes it monotone
    cum = np.cumsum([out[f"k{k}"] for k in range(1, KM + 1)])
    return dict(per_k=out, cumulative=[float(c) for c in cum])


def crossing(gaps_by_rung, per_list):
    """All interpolated total-n crossings of a gap sequence."""
    ns = [3 * p for p in per_list]            # total study n at 1:1:1
    g = [gaps_by_rung[p] for p in per_list]
    out = []
    for i in range(len(g) - 1):
        a, b = g[i], g[i + 1]
        if (a > 0 >= b) or (a < 0 <= b):
            t = a / (a - b)
            out.append(float(ns[i] + t * (ns[i + 1] - ns[i])))
    return out


def main(chrom=22):
    res = {"chrom": chrom, "n_splits": N_SPLITS, "share": {}, "crossing": {}}
    # fixed-composition 3-super-pop ladder (label = the per-sample SUPER-POP,
    # matching s51's fixed 1:1:1 composition ladder)
    G3, _, sup3, _ = s51._sup_label(SUPERPOP, chrom, 1500, 20000)
    per_seed_gap = {p: [] for p in LADDER}
    for p in LADDER:
        shares = []
        for s in range(N_SPLITS):
            r = share_ladder(G3, sup3, p, s)
            if r:
                shares.append(r["cumulative"])
            row = s48.cohort_cell(G3, sup3, per=p, seed=s)
            if row:
                per_seed_gap[p].append(row["gap"])
        if shares:
            arr = np.array(shares)
            res["share"][f"per{p}"] = dict(mean=arr.mean(0).tolist(),
                                           std=arr.std(0).tolist())
        print("per%-4d gap mean %+.3f sd %.3f (n=%d)" % (
            p, np.mean(per_seed_gap[p]), np.std(per_seed_gap[p]),
            len(per_seed_gap[p])), flush=True)
    mean_gap = {p: float(np.mean(per_seed_gap[p])) for p in LADDER}
    res["crossing"]["mean_sequence"] = crossing(mean_gap, LADDER)
    res["crossing"]["mean"] = (res["crossing"]["mean_sequence"][0]
                               if len(res["crossing"]["mean_sequence"]) == 1 else None)
    per_seed_cross = {}
    for s in range(N_SPLITS):
        gs = {}
        ok = True
        for p in LADDER:
            g = s48.cohort_cell(G3, sup3, per=p, seed=s)
            if not g:
                ok = False; break
            gs[p] = g["gap"]
        if ok:
            per_seed_cross[s] = crossing(gs, LADDER)
    res["crossing"]["per_seed_by_id"] = per_seed_cross
    single = [c[0] for c in per_seed_cross.values() if len(c) == 1]
    res["crossing"]["n_single_crossing_seeds"] = len(single)
    res["crossing"]["mean"] = float(np.mean(single)) if single else None
    res["crossing"]["sd"] = float(np.std(single)) if single else None
    res["crossing"]["mean_gap_by_rung"] = mean_gap
    with open(os.path.join(OUT, "s60_share_crossing.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("crossing mean", res["crossing"]["mean"],
          "per-seed", res["crossing"]["per_seed_by_id"],
          "sd", res["crossing"]["sd"], flush=True)
    return res


if __name__ == "__main__":
    t0 = time.time()
    payload = main()
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "s60_share_crossing.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
