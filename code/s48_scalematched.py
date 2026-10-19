"""S48 wrap-up - small-n structure ladder at matched PER-GROUP size.

The s48 ladder confounds structure with study size (single-pop n=28-36, 2-pop
~52-66, 3/5-pop 450/715). This reruns the same `cohort_cell` at a fixed
per-group sample n=20 for every rung, so the only thing that changes across the
ladder is the number of groups (i.e. structure); the TOTAL study n still grows
with the group count (single-pop n=20, 3-pop n=60, 5-pop n=100), which is
recorded per rung so the residual size confound is explicit.

Reports mean and std over the split seeds, plus the total study n, for the AF
channel, the oriented PC channel, the gap, and the two basis controls.

Output: results/raw/s48_scalematched.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
from real_multi import SUPERPOP, SUPERPOP5


OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
PER = 20
N_SPLITS = 5


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


def run(chrom=22, ld_window=1500, max_snps=20000):
    samples, sups, pops = s48.pop_table()
    G, lab, pos = s48.load_for(SUPERPOP5, chrom, ld_window, max_snps)
    order = [s for s, k in zip(samples, sups) if k in set(SUPERPOP5)]
    order_sup = np.array([sups[samples.index(s)] for s in order])
    res = {}
    for pop in ("GBR", "CHB"):
        m = lab == pop
        rows = [s48.cohort_cell(G[m], lab[m], per=PER, seed=s) for s in range(N_SPLITS)]
        res[f"single:{pop}"] = _agg([r for r in rows if r])
    for name, sp_keep, labarr in (("3pop", SUPERPOP, order_sup),
                                  ("5pop", SUPERPOP5, order_sup)):
        m = np.isin(labarr, sp_keep)
        rows = [s48.cohort_cell(G[m], labarr[m], per=PER, seed=s) for s in range(N_SPLITS)]
        res[name] = _agg([r for r in rows if r])
    payload = dict(chrom=chrom, per_group=PER, n_splits=N_SPLITS, cohorts=res)
    with open(os.path.join(OUT, "s48_scalematched.json"), "w") as f:
        json.dump(payload, f, indent=2)
    for name, a in res.items():
        print("%-12s ngrp %d n_tot %d  fst %.4f  AF %.3f  PC %.3f  gap %+.3f (sd %.3f)  "
              "refb %.3f  eigo %.3f" % (
                  name, a["n_group"], a["n_study"], a["fst_like"], a["af_auc"],
                  a["pc_auc"], a["gap"], a["gap_std"], a["pc_refbasis"], a["pc_eigonly"]),
              flush=True)
    return payload


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    run(a.chrom)
