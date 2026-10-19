"""S0b - sanity aggregation over seeds, for the lambda_GC figure in the paper.
Reports median-over-seeds uncorrected/corrected lambda_GC and attack sanity.
Output: results/raw/s0b_sanity_seeds.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    simulate_population, split_population3, release_objects, marginal_association,
    genomic_control_lambda, allele_freq, af_lrt_attack, anc_only_attack,
    auroc, pca,
)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def run(seeds=range(8)):
    unc, cor, atk, anc, noise = [], [], [], [], []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, nonmem = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        mask = rel["mask"]
        _, pu = marginal_association(site["G"], site["case"], 0, False)
        _, pc = marginal_association(site["G"], site["case"], 10, True)
        unc.append(genomic_control_lambda(pu))
        cor.append(genomic_control_lambda(pc))
        af_s = allele_freq(site["G"]); af_r = allele_freq(ref["G"])
        sm = af_lrt_attack(site["G"], af_s, af_r, mask=mask)
        sn = af_lrt_attack(nonmem["G"], af_s, af_r, mask=mask)
        lab = np.concatenate([np.ones(len(sm)), np.zeros(len(sn))])
        atk.append(auroc(np.concatenate([sm, sn]), lab))
        _, _, load_s = pca(site["G"], k=10)
        mean_s = site["G"][:, mask].mean(0)
        am = anc_only_attack(site["G"], load_s, mean_s, mask)
        an = anc_only_attack(nonmem["G"], load_s, mean_s, mask)
        anc.append(auroc(np.concatenate([am, an]), lab))
        rng = np.random.default_rng(seed + 7)
        gm = rng.integers(0, 3, size=site["G"].shape).astype(np.int8)
        gn = rng.integers(0, 3, size=nonmem["G"].shape).astype(np.int8)
        xm = af_lrt_attack(gm, af_s, af_r, mask=mask)
        xn = af_lrt_attack(gn, af_s, af_r, mask=mask)
        noise.append(auroc(np.concatenate([xm, xn]), lab))
    out = dict(
        seeds=list(seeds),
        lambda_uncorrected_mean=float(np.mean(unc)),
        lambda_uncorrected_median=float(np.median(unc)),
        lambda_uncorrected_range=[float(np.min(unc)), float(np.max(unc))],
        lambda_corrected_mean=float(np.mean(cor)),
        lambda_corrected_median=float(np.median(cor)),
        attack_auroc_mean=float(np.mean(atk)),
        study_pc_distance_auroc_mean=float(np.mean(anc)),
        noise_auroc_mean=float(np.mean(noise)),
    )
    print(json.dumps(out, indent=2))
    with open(os.path.join(OUT, "s0b_sanity_seeds.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out


if __name__ == "__main__":
    run()
