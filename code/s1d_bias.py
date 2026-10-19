"""S1d - Quantify the self-reference bias reproducibly.

On the SAME panels, compare the biased protocol (attacker's own reference
individuals used as non-members) against the disjoint protocol (independent
non-member panel), at F5_full. Reports the AUROC inflation.

Output: results/raw/s1d_bias.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    simulate_population, split_population3, release_objects, allele_freq,
    af_lrt_attack, auroc,
)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def run(seeds=(0, 1, 2, 3, 4)):
    rows = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, nonmem = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        mask = rel["mask"]
        af_s = allele_freq(site["G"]); af_r = allele_freq(ref["G"])
        sm = af_lrt_attack(site["G"], af_s, af_r, mask=mask)
        # biased: negatives = attacker reference (pulled toward ref_af)
        sref = af_lrt_attack(ref["G"], af_s, af_r, mask=mask)
        lab_b = np.concatenate([np.ones(len(sm)), np.zeros(len(sref))])
        a_bias = auroc(np.concatenate([sm, sref]), lab_b)
        # disjoint: negatives = independent non-member panel
        snm = af_lrt_attack(nonmem["G"], af_s, af_r, mask=mask)
        lab_d = np.concatenate([np.ones(len(sm)), np.zeros(len(snm))])
        a_disc = auroc(np.concatenate([sm, snm]), lab_d)
        rows.append(dict(seed=seed, biased=float(a_bias), disjoint=float(a_disc),
                         inflation=float(a_bias - a_disc)))
        print(rows[-1], flush=True)
    out = dict(rows=rows,
               biased_mean=float(np.mean([r["biased"] for r in rows])),
               disjoint_mean=float(np.mean([r["disjoint"] for r in rows])),
               inflation_mean=float(np.mean([r["inflation"] for r in rows])))
    print(json.dumps(out, indent=2))
    with open(os.path.join(OUT, "s1d_bias.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out


if __name__ == "__main__":
    run()
