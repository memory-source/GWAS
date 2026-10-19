"""S55 - real within-population ordering across marker density.

s30 shows the PC-vs-AF ordering on the real MIXED (3-super-population) cohort is
density-stable across a 50x thinning range. The real within-population ladder
(s48/s50) runs at a single density (15 kb window). This runs the s48 within-
population cohorts (GBR/CHB/YRI, chr22) at two marker densities so the
within-population PC>AF ordering is shown density-stable too, closing the
"one-thinning-contrast" caveat on the real-data M axis.

Densities: ld_window 1500 bp (~3.4k markers) and 300 bp (~22k markers, near all
common SNPs). Attack code identical to s48 (`cohort_cell`).

Output: results/raw/s55_real_density.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
from real_multi import SUPERPOP5

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
DENSITIES = {"D_thin": 1500, "D_med": 3000, "D_dense": 300}
N_SPLITS = 5


def run(chrom=22, pops=("GBR", "CHB", "YRI")):
    res = {"chrom": chrom, "n_splits": N_SPLITS, "cells": {}}
    for name, lw in DENSITIES.items():
        G, lab, _ = s48.load_for(SUPERPOP5, chrom, ld_window=lw, max_snps=500000)
        for pop in pops:
            m = lab == pop
            if m.sum() < 12:
                continue
            rows = [s48.cohort_cell(G[m], lab[m], per=150, seed=s) for s in range(N_SPLITS)]
            rows = [r for r in rows if r]
            if not rows:
                continue
            agg = {}
            for k in rows[0]:
                if not isinstance(rows[0][k], (int, float)):
                    continue
                vals = np.array([r[k] for r in rows], float)
                agg[k] = float(vals.mean())
                if k not in ("n_study", "n_group", "M"):
                    agg[k + "_std"] = float(vals.std())
            res["cells"][f"{name}:{pop}"] = agg
            print("%-10s %-4s M %5d  AF %.3f  PC %.3f  gap %+.3f (sd %.3f)  refb %.3f  eigo %.3f" % (
                name, pop, agg["M"], agg["af_auc"], agg["pc_auc"], agg["gap"],
                agg["gap_std"], agg["pc_refbasis"], agg["pc_eigonly"]), flush=True)
    with open(os.path.join(OUT, "s55_real_density.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    t0 = time.time()
    payload = run()
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "s55_real_density.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
