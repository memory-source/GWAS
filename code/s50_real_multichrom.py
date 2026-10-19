"""S50 - real structure ladder on independent chromosomes.

Replicates s48's structure ladder (within-population PC>AF -> across-population
AF>PC) on chr1 and chr2, using the cached panels. This bounds the real-data law
to more than one chromosome.

Output: results/raw/s50_real_structure_multichrom.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")

KEEP = ["pop:GBR", "pop:CHB", "pop:YRI", "super:EUR-2pop",
        "super:3pop", "super:5pop"]


def run(chroms=(1, 2), n_splits=5):
    out = {}
    for chrom in chroms:
        res = s48.run_chrom(chrom, n_splits=n_splits)
        out[f"chr{chrom}"] = {k: res[k] for k in res if k in KEEP}
        for name in sorted(out[f"chr{chrom}"], key=lambda k: out[f"chr{chrom}"][k]["fst_like"]):
            a = out[f"chr{chrom}"][name]
            print("chr%d %-16s fst %.4f AF %.3f PC %.3f gap %+.3f refb %.3f eigo %.3f" % (
                chrom, name, a["fst_like"], a["af_auc"], a["pc_auc"], a["gap"],
                a["pc_refbasis"], a["pc_eigonly"]), flush=True)
        with open(os.path.join(OUT, "s50_real_structure_multichrom.json"), "w") as f:
            json.dump(dict(per_chrom=out, n_splits=n_splits), f, indent=2)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default="1,2")
    ap.add_argument("--splits", type=int, default=5)
    a = ap.parse_args()
    run(tuple(int(c) for c in a.chroms.split(",")), a.splits)
