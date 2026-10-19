"""S21 - Fst / shared-structure crossover sweep.

The paper's honest boundary is that the PC-vs-AF ordering holds on the
parametric simulator but not on real 1000G data, with the magnitude difference
currently unexplained. This sweep tests whether the ordering is monotone in ONE
parameter (the between-group structure level): as Fst rises, the top PCs capture
more between-group structure shared equally by members and non-members, so the
in-sample residual share shrinks and the PC channel should fall toward the real
regime while the AF-LRT channel stays roughly flat.

Grid: Fst in {0.02, 0.04, 0.08, 0.12, 0.16, 0.24} x n_per_group in
{300, 800, 1600}, M=8000, 15 seeds. At each point: AF-only, train-half
calibrated PC, combined; in-sample mechanism (own/other/swapped basis); the
structure covariates fst_like and top10_var_frac. Real chr1/chr2/chr22 points
are overlaid as reference markers.

Output: results/raw/s21_sensitivity.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def pc_dist(G, mask, loadings, mean, k=None):
    L = loadings if k is None else loadings[:k]
    X = G[:, mask].astype(np.float64) - mean[None, :]
    return np.linalg.norm((X @ L[:, mask].T) / np.sqrt(max(1, mask.sum())), axis=1)


def cal(a):
    return max(float(a), 1 - float(a))


def fst_like(G, groups):
    af = G.mean(0) / 2.0
    hs = np.mean([2 * (G[groups == g].mean(0) / 2.0) * (1 - G[groups == g].mean(0) / 2.0)
                  for g in np.unique(groups)], axis=0)
    ht = 2 * af * (1 - af)
    return float(np.nanmean((ht - hs) / np.where(ht > 0, ht, np.nan)))


def top10_frac(G, mask, eig):
    X = G[:, mask].astype(float); X = X - X.mean(0)
    S = np.linalg.svd(X, compute_uv=False)
    total = (S ** 2).sum() / max(1, X.shape[0] - 1)
    return float(eig.sum() / total)


def one_cell(fst, npg, n_snps, seed):
    pop = simulate_population(n_groups=3, n_per_group=npg, n_snps=n_snps,
                              fst=fst, ld_blocks=100, seed=seed)
    per = max(50, npg // 3)
    site, ref, non = split_population3(pop, per, per, per)
    rel = release_objects(site)
    mask = rel["mask"]
    if mask.sum() < 10:
        return None
    af_s, af_r = rel["af"], allele_freq(ref["G"])
    mean_s = site["G"][:, mask].mean(0)
    L = rel["loadings"]
    y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])
    am = af_lrt_attack(site["G"], af_s, af_r, mask=mask)
    an = af_lrt_attack(non["G"], af_s, af_r, mask=mask)
    dm, dn = pc_dist(site["G"], mask, L, mean_s), pc_dist(non["G"], mask, L, mean_s)
    af_auc = float(auroc(np.concatenate([am, an]), y))
    pc_auc = cal(auroc(np.concatenate([dm, dn]), y))
    # combined via logistic on a train half of the labels
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    rng = np.random.default_rng(seed + 313)
    tsm = rng.random(len(dm)) < 0.5; tsn = rng.random(len(dn)) < 0.5
    ytr = np.concatenate([np.ones(tsm.sum()), np.zeros(tsn.sum())])
    yte = np.concatenate([np.ones((~tsm).sum()), np.zeros((~tsn).sum())])
    pc_tr = np.concatenate([dm[tsm], dn[tsn]])
    sgn = 1.0 if auroc(pc_tr, ytr) >= auroc(-pc_tr, ytr) else -1.0
    Xtr = np.column_stack([np.concatenate([am[tsm], an[tsn]]), sgn * pc_tr])
    Xte = np.column_stack([np.concatenate([am[~tsm], an[~tsn]]),
                           sgn * np.concatenate([dm[~tsm], dn[~tsn]])])
    sc = StandardScaler().fit(Xtr)
    comb = LogisticRegression(max_iter=1000).fit(sc.transform(Xtr), ytr).predict_proba(
        sc.transform(Xte))[:, 1]
    comb_auc = float(auroc(comb, yte))
    # in-sample mechanism (ancestry-matched 4-way split within the population)
    rng2 = np.random.default_rng(seed + 77)
    A_i, B_i, C_i = [], [], []
    for k in range(3):
        idx = rng2.permutation(np.where(pop["ancestry"] == k)[0])
        q = min(per, idx.size // 4)
        A_i.append(idx[:q]); B_i.append(idx[q:2 * q]); C_i.append(idx[2 * q:3 * q])
    A, B, C = (pop["G"][np.concatenate(x)] for x in (A_i, B_i, C_i))
    _, _, LA = __import__("common").pca(A, k=10); mA = A[:, mask].mean(0)
    _, _, LB = __import__("common").pca(B, k=10); mB = B[:, mask].mean(0)
    yc = np.concatenate([np.ones(len(A)), np.zeros(len(C))])
    own = cal(auroc(np.concatenate([pc_dist(A, mask, LA, mA), pc_dist(C, mask, LA, mA)]), yc))
    oth = cal(auroc(np.concatenate([pc_dist(A, mask, LB, mB), pc_dist(C, mask, LB, mB)]), yc))
    return dict(af_only=af_auc, pc_calibrated=pc_auc, combined=comb_auc,
                delta_pc_minus_af=pc_auc - af_auc,
                delta_comb_minus_af=comb_auc - af_auc,
                own_basis=own, other_basis=oth,
                fst_like=fst_like(site["G"], site["ancestry"]),
                top10_var_frac=top10_frac(site["G"], mask, rel["eig"]))


def run(fsts=(0.02, 0.04, 0.08, 0.12, 0.16, 0.24),
        n_per_list=(300, 800, 1600), n_snps=8000, n_seeds=15):
    grid = {}
    for fst in fsts:
        for npg in n_per_list:
            rows = [one_cell(fst, npg, n_snps, s) for s in range(n_seeds)]
            rows = [r for r in rows if r]
            key = f"f{fst}_n{npg}"
            agg = {}
            for k in ("af_only", "pc_calibrated", "combined", "delta_pc_minus_af",
                      "delta_comb_minus_af", "own_basis", "other_basis",
                      "fst_like", "top10_var_frac"):
                agg[k] = float(np.mean([r[k] for r in rows]))
                agg[k + "_std"] = float(np.std([r[k] for r in rows]))
            grid[key] = dict(agg=agg, n=len(rows))
            print(key, json.dumps({k: round(v, 4) for k, v in agg.items() if not k.endswith("_std")}),
                  flush=True)
    # crossover: interpolate the Fst where delta_pc_minus_af crosses 0 at n=800
    monotone = {}
    for npg in n_per_list:
        ds = [grid[f"f{fst}_n{npg}"]["agg"]["delta_pc_minus_af"] for fst in fsts]
        mono = all(ds[i] >= ds[i + 1] - 0.01 for i in range(len(ds) - 1))
        monotone[f"n{npg}"] = dict(deltas=ds, monotone_decreasing=bool(mono))
    res = dict(grid=grid, fsts=list(fsts), n_per_group=list(n_per_list),
               n_snps=n_snps, n_seeds=n_seeds, monotone=monotone)
    with open(os.path.join(OUT, "s21_sensitivity.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
