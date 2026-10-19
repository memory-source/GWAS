"""Paired negation sweep for the two SLOW cells only (n=1600/M=8000 and
n=400/M=2000/Fst=0.24), 200 seeds each, writing per-seed deltas so the summary
can pool them with the fast cells. Run as ONE process (the slow cell is
BLAS-heavy; parallel shards thrash a shared box)."""
import json
import os
import sys
import time
import numpy as np

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
import s29_negation as S

OUT = "/data/lvvyongji1/GWAS/results/raw"
SLOW = S.PAIRED_SLOW

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    t0 = time.time()
    res = {"paired": S.paired_sweep(SLOW, n_seeds=n, seed0=1000)}
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "s29_negation_slow.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("slow paired done", res["elapsed_s"], "s", flush=True)
