"""Re-run the s29 ladders only (sim/real/msprime exact k=1..10 + eig-only),
writing the ladder blocks into a dedicated shard file. Kept separate so we can
extend the ladder to larger k without re-running the 600-seed paired sweep."""
import json
import os
import sys
import time

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
import s29_negation as S

OUT = "/data/lvvyongji1/GWAS/results/raw"

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    res = {}
    if which in ("all", "sim"):
        res["ladder_sim"] = S.ladder_sim()
    if which in ("all", "real"):
        res["ladder_real"] = S.ladder_real()
    if which in ("all", "msprime"):
        res["ladder_msprime"] = S.ladder_msprime()
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "s29_ladder.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("ladders written", which, res["elapsed_s"], "s", flush=True)
