"""S51 - real low-structure / large-n ladder (round-15 limitation closure).

Closes the paper's "no real low-structure large-n cell" limitation. Two arms on
real 1000G phase-3:

  Arm A (pooled within-super-population): one super-population (EUR/EAS/AFR),
  the whole pool as a single group, disjoint 3-way study/reference/non split at
  the largest per-group n the pool allows. Measured between-group structure = 0
  by construction; this is the largest real low-structure cell 1000G supports.

  Arm B (fixed-composition size ladder): the 3-super-population cohort with the
  composition HELD FIXED (AFR+EAS+EUR) and the per-group n swept
  {40,80,120,160}, so only the total study size changes across rungs. This
  isolates the total-n dependence of the PC-vs-AF gap from composition.

Per cell: af_auc, pc_auc (held-out-half oriented), gap, measured fst_like, the
reference-fitted-basis and eigenvalue-only controls. 5 split seeds each.

Output: results/raw/s51_real_lowstruct.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
from real_multi import SUPERPOP, SUPERPOP5


OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
N_SPLITS = 5
POOL_PER = 100000          # effectively unbounded: cohort_cell caps at pool//3
LADDER = (40, 80, 120, 160)


def _agg(rows):
    out = {}
    for k in rows[0]:
        if not isinstance(rows[0][k], (int, float)):
            continue
        vals = np.array([r[k] for r in rows], float)
        out[k] = float(vals.mean())
        if k not in ("n_study", "n_group", "M"):
            out[k + "_std"] = float(vals.std())
    return out


def _sup_label(sups_keep, chrom, ld_window=1500, max_snps=20000):
    G, lab, pos = s48.load_for(sups_keep, chrom, ld_window, max_snps)
    samples, sups, pops = s48.pop_table()
    order = [s for s, k in zip(samples, sups) if k in set(sups_keep)]
    sup_arr = np.array([sups[samples.index(s)] for s in order])
    return G, lab, sup_arr, pos


def run(chrom=22, ld_window=1500, max_snps=20000):
    res = {"chrom": chrom, "n_splits": N_SPLITS, "pooled": {}, "ladder": {}}
    # Arm A: pooled within-super-population (5-super-pop load so every pop exists)
    G5, _, sup5, _ = _sup_label(SUPERPOP5, chrom, ld_window, max_snps)
    for sp in ("EUR", "EAS", "AFR"):
        m = sup5 == sp
        n_pool = int(m.sum())
        rows = [s48.cohort_cell(G5[m], np.array(["pool"] * n_pool), per=POOL_PER, seed=s)
                for s in range(N_SPLITS)]
        rows = [r for r in rows if r]
        if rows:
            a = _agg(rows)
            a["n_pool"] = n_pool
            res["pooled"][sp] = a
            print("pooled %-4s n_pool %d  n_study %d  AF %.3f  PC %.3f  gap %+.3f "
                  "(sd %.3f)  refb %.3f  eigo %.3f" % (
                      sp, n_pool, a["n_study"], a["af_auc"], a["pc_auc"], a["gap"],
                      a["gap_std"], a["pc_refbasis"], a["pc_eigonly"]), flush=True)
    # Arm B: fixed-composition 3-super-pop size ladder
    G3, _, sup3, _ = _sup_label(SUPERPOP, chrom, ld_window, max_snps)
    for per in LADDER:
        rows = [s48.cohort_cell(G3, sup3, per=per, seed=s) for s in range(N_SPLITS)]
        rows = [r for r in rows if r]
        if rows:
            a = _agg(rows)
            res["ladder"][f"per{per}"] = a
            print("ladder per %3d  n_tot %3d  fst %.4f  AF %.3f  PC %.3f  gap %+.3f "
                  "(sd %.3f)  refb %.3f  eigo %.3f" % (
                      per, a["n_study"], a["fst_like"], a["af_auc"], a["pc_auc"],
                      a["gap"], a["gap_std"], a["pc_refbasis"], a["pc_eigonly"]),
                  flush=True)
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s51_real_lowstruct_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
