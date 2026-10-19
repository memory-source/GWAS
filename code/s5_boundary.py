"""S5 - Population-level boundary curve (the actionable deliverable).

For each site size and release precision, measure the *decision threshold* at
which an attacker can call membership: here, per-individual attack scores are
aggregated into a study-level AUC, and we report the population-level decision
boundary rather than a single point.

Fast (2 seeds, differential_evolution-free: coarse grid). Output:
results/raw/s5_boundary.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    simulate_population, split_population3, release_objects, allele_freq,
    af_lrt_attack, auroc, tpr_at_fpr,
)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")
os.makedirs(OUT, exist_ok=True)


def run(seeds=(0, 1)):
    grid = {}
    for n_site in (50, 100, 200, 400, 800):
        for dec in (1, 2, None):
            aucs, tprs = [], []
            for seed in seeds:
                pop = simulate_population(n_groups=3, n_per_group=n_site + 600,
                                          n_snps=4000, seed=seed)
                site, ref, nonmem = split_population3(pop, n_site, 300, 300)
                rel = release_objects(site)
                af = rel["af"] if dec is None else np.round(rel["af"], dec)
                r_af = allele_freq(ref["G"])
                sm = af_lrt_attack(site["G"], af, r_af, mask=rel["mask"])
                sn = af_lrt_attack(nonmem["G"], af, r_af, mask=rel["mask"])
                lab = np.concatenate([np.ones(len(sm)), np.zeros(len(sn))])
                sc = np.concatenate([sm, sn])
                aucs.append(auroc(sc, lab))
                tprs.append(tpr_at_fpr(sc, lab, 0.01))
            key = f"n{n_site}_dec{dec}"
            grid[key] = dict(n_site=n_site, decimals=dec,
                             auroc=float(np.mean(aucs)), tpr_at_fpr1=float(np.mean(tprs)))
            print(key, json.dumps(grid[key]), flush=True)
    with open(os.path.join(OUT, "s5_boundary.json"), "w") as f:
        json.dump(grid, f, indent=2)
    return grid


if __name__ == "__main__":
    run()
