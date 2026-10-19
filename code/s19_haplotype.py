"""S19 - haplotype/LD channel: a genuinely distinct attack family.

All other channels are marginal-frequency attacks (AF-LRT, correlation,
learned deviation all use per-SNP marginals). This one uses the *joint*
distribution of adjacent marker pairs — the within-window LD/haplotype
structure — which the marginal AF-LRT cannot see.

Score: for each adjacent marker pair (j, j+1), estimate the joint genotype
distribution P(g_j, g_{j+1}) from the study panel and from the reference panel
(a 3x3 table, smoothed), then score a target by
  sum_j [ log P_study(g_j, g_{j+1}) - log P_ref(g_j, g_{j+1}) ].
Members (whose genotypes shaped the study's joint LD table) score higher.

Output: results/raw/s19_haplotype.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def _pair_tables(G):
    """Joint genotype tables for adjacent marker pairs: list of 3x3 counts."""
    n, m = G.shape
    tabs = []
    for j in range(m - 1):
        t = np.zeros((3, 3))
        np.add.at(t, (G[:, j], G[:, j + 1]), 1)
        tabs.append(t)
    return tabs


def _loglik(tabs, G):
    """Per-individual log-likelihood under a set of adjacent-pair tables."""
    ll = np.zeros(G.shape[0])
    for j, t in enumerate(tabs):
        p = (t + 0.5) / (t.sum() + 4.5)
        ll += np.log(p[G[:, j], G[:, j + 1]] + 1e-12)
    return ll


def hap_ld_attack(target_G, study_G, ref_G, mask):
    idx = np.where(mask)[0]
    st = _pair_tables(study_G[:, idx])
    rt = _pair_tables(ref_G[:, idx])
    return _loglik(st, target_G[:, idx]) - _loglik(rt, target_G[:, idx])


def run(seeds=(0, 1, 2, 3, 4)):
    rows = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, non = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        mask = rel["mask"]
        y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])
        hm = hap_ld_attack(site["G"], site["G"], ref["G"], mask)
        hn = hap_ld_attack(non["G"], site["G"], ref["G"], mask)
        hap_auc = float(auroc(np.concatenate([hm, hn]), y))
        af_s, af_r = rel["af"], allele_freq(ref["G"])
        am = af_lrt_attack(site["G"], af_s, af_r, mask=mask)
        an = af_lrt_attack(non["G"], af_s, af_r, mask=mask)
        af_auc = float(auroc(np.concatenate([am, an]), y))
        # Out-of-sample transfer: fit BOTH joint tables on panels that saw
        # NEITHER the study nor the targets (the attacker reference panel and a
        # fresh independent population D2), then score the same members vs
        # non-members. If the channel is in-sample-ness of the released table,
        # this is at chance.
        pop2 = simulate_population(n_groups=3, n_per_group=300, n_snps=4000, seed=seed + 5000)
        D2 = pop2["G"]
        hm_oos = hap_ld_attack(site["G"], ref["G"], D2, mask)
        hn_oos = hap_ld_attack(non["G"], ref["G"], D2, mask)
        oos = float(auroc(np.concatenate([hm_oos, hn_oos]), y))
        row = dict(seed=seed, hap_ld=hap_auc, af_only=af_auc, hap_oos_transfer=oos)
        rows.append(row)
        print(json.dumps({k: round(v, 4) if isinstance(v, float) else v for k, v in row.items()}),
              flush=True)
    keys = ["hap_ld", "af_only", "hap_oos_transfer"]
    agg = {k: (float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows])))
           for k in keys}
    res = dict(per_seed=rows, agg=agg)
    with open(os.path.join(OUT, "s19_haplotype.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("agg", json.dumps({k: round(v[0], 4) for k, v in agg.items()}), flush=True)
    return res


if __name__ == "__main__":
    run()
