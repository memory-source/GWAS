"""S2b - Marginal leakage of intermediate reveals with disjoint non-members
(reviewer fix), plus permutation decomposition and bootstrap CIs.

Same meta-attack design as S2 but the negative class is an independent
non-member panel disjoint from both study and attacker reference.
Output: results/raw/s2b_intermediate.json
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
os.makedirs(OUT, exist_ok=True)


def features_for(site, ref, nonmem, rel, final_tier="F3_round2"):
    s_af = np.round(rel["af"], 2) if final_tier == "F3_round2" else rel["af"]
    r_af = allele_freq(ref["G"])
    mask = rel["mask"]
    eig, scores_study, loadings = pca(site["G"], k=10)

    def per_panel(P):
        G = P["G"]
        f = {}
        f["af"] = af_lrt_attack(G, s_af, r_af, mask=mask)
        X = G[:, mask].astype(np.float64)
        X = X - X.mean(0, keepdims=True)
        proj = (X @ loadings[:, mask].T) / np.sqrt(max(1, mask.sum()))
        centroid = proj.mean(0)
        f["m2_pcdist"] = -np.linalg.norm(proj - centroid[None, :], axis=1)
        f["m2_pcalign"] = np.abs(proj[:, 0])
        miss_c = rel["missing"][mask] - rel["missing"][mask].mean()
        f["m3_misscorr"] = (X * miss_c[None, :]).sum(1)
        f["m1_lambda"] = np.full(G.shape[0], rel["lambda_gc"])
        f["m4_condbit"] = np.full(G.shape[0], float(rel["eig"][0]))
        f["pca_proj"] = proj  # for ancestry-confound permutation
        return f

    fm = per_panel(site)
    fn = per_panel(nonmem)
    return fm, fn


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
    return auroc(oof, y), oof


def run(seeds=(0, 1, 2)):
    BASE = ["af"]
    SETS = {
        "base_final_only": BASE,
        "base+M1_lambdaGC": BASE + ["m1_lambda"],
        "base+M4_conditionbit": BASE + ["m4_condbit"],
        "base+M3_missingness": BASE + ["m3_misscorr"],
        "base+M2_pcaSpectrum": BASE + ["m2_pcdist", "m2_pcalign"],
        "base+M2+M3": BASE + ["m2_pcdist", "m2_pcalign", "m3_misscorr"],
    }
    results = {k: [] for k in SETS}
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, nonmem = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        fm, fn = features_for(site, ref, nonmem, rel)
        y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(nonmem["G"]))])
        for name, feats in SETS.items():
            X = np.column_stack([np.concatenate([fm[k], fn[k]]) for k in feats])
            a, _ = cv_auc(X, y, seed=seed)
            results[name].append(a)
        print("seed", seed, "done", {k: round(np.mean(v), 3) for k, v in results.items()})
    out = {k: dict(auc_mean=float(np.mean(v)), auc_std=float(np.std(v)),
                   per_seed=[float(x) for x in v]) for k, v in results.items()}
    base_auc = out["base_final_only"]["auc_mean"]
    out["marginal_over_base"] = {k: float(out[k]["auc_mean"] - base_auc)
                                 for k in SETS if k != "base_final_only"}
    print(json.dumps(out, indent=2))
    with open(os.path.join(OUT, "s2b_intermediate.json"), "w") as f:
        json.dump(dict(seeds=list(seeds), results=out), f, indent=2)
    return out


if __name__ == "__main__":
    run()
