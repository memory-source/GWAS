"""S6 - Dimensionality sweep for the PCA-spectrum leak (round-2 must-do).

The PC-distance leak's magnitude depends on the ratio m/n (SNPs per member).
Sweep n_snps and n_per_group with everything else fixed, reporting the
pcdist-only AUROC and the base+M2 delta. Output: results/raw/s6_dim.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    simulate_population, split_population3, release_objects, allele_freq,
    af_lrt_attack, auroc, pca,
)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def cv_auc(X, y, seed=0, folds=5):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr])
        clf = LogisticRegression(max_iter=2000)
        clf.fit(sc.transform(X[tr]), y[tr])
        oof[te] = clf.predict_proba(sc.transform(X[te]))[:, 1]
    return auroc(oof, y)


def m2(G, mask, loadings, center_by):
    X = G[:, mask].astype(np.float64) - center_by[None, :]
    proj = (X @ loadings[:, mask].T) / np.sqrt(max(1, mask.sum()))
    return -np.linalg.norm(proj, axis=1)


def run(seeds=(0, 1)):
    rows = []
    for n_snps in (1000, 3000, 6000):
        for n_site in (100, 300, 600):
            a_pcd, a_base, a_bm = [], [], []
            for seed in seeds:
                pop = simulate_population(n_groups=3, n_per_group=n_site + 600,
                                          n_snps=n_snps, seed=seed, ld_blocks=0)
                site, ref, nonmem = split_population3(pop, n_site, 300, 300)
                rel = release_objects(site)
                mask = rel["mask"]
                r_af = allele_freq(ref["G"])
                s_af = np.round(rel["af"], 2)
                y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(nonmem["G"]))])
                mean_s = site["G"][:, mask].mean(0)
                _, _, load = pca(site["G"], k=10)
                pm = m2(site["G"], mask, load, mean_s)
                pn = m2(nonmem["G"], mask, load, mean_s)
                bm = af_lrt_attack(site["G"], s_af, r_af, mask=mask)
                bn = af_lrt_attack(nonmem["G"], s_af, r_af, mask=mask)
                a_pcd.append(cv_auc(np.concatenate([pm, pn])[:, None], y, seed))
                a_base.append(cv_auc(np.concatenate([bm, bn])[:, None], y, seed))
                X = np.column_stack([np.concatenate([bm, bn]), np.concatenate([pm, pn])])
                a_bm.append(cv_auc(X, y, seed))
            row = dict(n_snps=n_snps, n_site=n_site, ratio=round(n_snps / n_site, 1),
                       pcdist_only=float(np.mean(a_pcd)),
                       base=float(np.mean(a_base)),
                       base_M2=float(np.mean(a_bm)),
                       delta=float(np.mean(a_bm) - np.mean(a_base)))
            rows.append(row); print(json.dumps(row), flush=True)
    with open(os.path.join(OUT, "s6_dim.json"), "w") as f:
        json.dump(dict(rows=rows), f, indent=2)
    return rows


if __name__ == "__main__":
    run()
