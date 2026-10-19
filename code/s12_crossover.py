"""S12 - PC-vs-AF crossover map over site size and marker count.

Result-to-claim gap: "dominant intermediate" and the +0.10 increment are
operating-point-specific (S6 showed +0.006..+0.136). Here we map, at each
(n_per_group, n_snps), the paired difference between the study-fitted PC
channel and the AF-only channel, with the PC orientation calibrated on a
held-out train half. This locates where the PC channel overtakes AF and where
the increment collapses.

Output: results/raw/s12_crossover.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def cell(n_per_group, n_snps, n_seeds=12):
    afv, pcv, combv = [], [], []
    daf, dpc, dcomb = [], [], []
    for seed in range(n_seeds):
        pop = simulate_population(n_groups=3, n_per_group=n_per_group,
                                  n_snps=n_snps, seed=seed)
        per = max(10, int(n_per_group * 0.3))
        site, ref, non = split_population3(pop, per, per, per)
        rel = release_objects(site)
        mask = rel["mask"]
        if mask.sum() < 5:
            continue
        af_s, af_r = rel["af"], allele_freq(ref["G"])
        mean_s = site["G"][:, mask].mean(0)
        load_s = rel["loadings"]

        def pc(G):
            X = G[:, mask].astype(np.float64) - mean_s[None, :]
            proj = (X @ load_s[:, mask].T) / np.sqrt(max(1, mask.sum()))
            return np.linalg.norm(proj, axis=1)

        dm, dn = pc(site["G"]), pc(non["G"])
        am = af_lrt_attack(site["G"], af_s, af_r, mask=mask)
        an = af_lrt_attack(non["G"], af_s, af_r, mask=mask)
        rng = np.random.default_rng(seed + 7717)
        tsm = rng.random(len(dm)) < 0.5
        tsn = rng.random(len(dn)) < 0.5
        ytr = np.concatenate([np.ones(tsm.sum()), np.zeros(tsn.sum())])
        yte = np.concatenate([np.ones((~tsm).sum()), np.zeros((~tsn).sum())])
        pc_tr = np.concatenate([dm[tsm], dn[tsn]])
        pc_te = np.concatenate([dm[~tsm], dn[~tsn]])
        sgn = 1.0 if auroc(pc_tr, ytr) >= auroc(-pc_tr, ytr) else -1.0
        af_te = np.concatenate([am[~tsm], an[~tsn]])
        a_af = float(auroc(af_te, yte))
        a_pc = float(auroc(sgn * pc_te, yte))
        afv.append(a_af); pcv.append(a_pc)
        daf.append(a_af); dpc.append(a_pc)
        # combined via logistic regression fit on train
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        sc = StandardScaler().fit(np.column_stack([np.concatenate([am[tsm], an[tsn]]),
                                                   sgn * pc_tr]))
        clf = LogisticRegression(max_iter=2000).fit(
            sc.transform(np.column_stack([np.concatenate([am[tsm], an[tsn]]), sgn * pc_tr])), ytr)
        comb = clf.predict_proba(sc.transform(np.column_stack([af_te, sgn * pc_te])))[:, 1]
        a_comb = float(auroc(comb, yte))
        combv.append(a_comb); dcomb.append(a_comb - a_af)
    return dict(af_only=float(np.mean(afv)), pc_calibrated=float(np.mean(pcv)),
                combined=float(np.mean(combv)),
                delta_pc_minus_af=float(np.mean(np.array(pcv) - np.array(afv))),
                delta_pc_minus_af_std=float(np.std(np.array(pcv) - np.array(afv))),
                delta_comb_minus_af=float(np.mean(np.array(combv) - np.array(afv))),
                delta_comb_minus_af_std=float(np.std(np.array(combv) - np.array(afv))),
                n_kept=len(afv))


def run():
    n_per_list = [400, 800, 1600]
    n_snp_list = [500, 1000, 2000, 4000, 8000]
    grid = {}
    for npg in n_per_list:
        for nsn in n_snp_list:
            c = cell(npg, nsn)
            grid[f"n{npg}_m{nsn}"] = c
            print(f"n={npg} m={nsn}", json.dumps({k: round(v, 4) for k, v in c.items()}),
                  flush=True)
    res = dict(grid=grid, n_per_group=n_per_list, n_snps=n_snp_list, sim_seeds=12)
    with open(os.path.join(OUT, "s12_crossover.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
