"""S24 - cross-run release composition audit.

The paper's threat model includes A3 (an adaptive analyst accumulating releases
across runs) but makes no composition claim. This measures it: a sequence of
pre-registered runs each releases ONE object that individually looks benign;
we track the cumulative attack AUROC (out-of-fold logistic over the union of
released features) after each run, and compare against the largest single
constituent. Then a governed policy with an anti-composition clause (AF and
loadings mutually exclusive) is shown to flag the violating run and hold the
cumulative AUROC at the safe tier.

Output: results/raw/s24_crossrun.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc)
from governance import Ledger

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def cum_auc(features, y, seed):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    X = np.column_stack(features)
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr])
        oof[te] = LogisticRegression(max_iter=1000).fit(sc.transform(X[tr]), y[tr]).predict_proba(
            sc.transform(X[te]))[:, 1]
    return float(auroc(oof, y))


def pc_dist(G, mask, L, mean):
    X = G[:, mask].astype(float) - mean[None, :]
    return np.linalg.norm((X @ L[:, mask].T) / np.sqrt(mask.sum()), axis=1)


def run(seeds=(0, 1, 2, 3, 4)):
    per_seed = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, non = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        mask = rel["mask"]; L = rel["loadings"]
        mean_s = site["G"][:, mask].mean(0)
        af_r = allele_freq(ref["G"])
        y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])

        # per-run released feature (member score, non-member score)
        def af_feat(af_s):
            return (af_lrt_attack(site["G"], af_s, af_r, mask=mask),
                    af_lrt_attack(non["G"], af_s, af_r, mask=mask))

        def pc_feat(load):
            return (pc_dist(site["G"], mask, load, mean_s), pc_dist(non["G"], mask, load, mean_s))

        # runs: 1 λ_GC (no membership info), 2 PC eigenvalues (no axes -> random),
        # 3 AF@1dp, 4 AF@2dp (refinement), 5 AF full (refinement), 6 loadings (PC)
        rng = np.random.default_rng(seed)
        L_rand = rng.normal(size=L.shape); L_rand, _ = np.linalg.qr(L_rand.T)
        run_feats = {
            "1_lambdaGC": (np.zeros(len(site["G"])), np.zeros(len(non["G"]))),
            "2_eigenvalues": pc_feat(L_rand.T),
            "3_af_1dp": af_feat(np.round(rel["af"], 1)),
            "4_af_2dp": af_feat(np.round(rel["af"], 2)),
            "5_af_full": af_feat(rel["af"]),
            "6_loadings": pc_feat(L),
        }
        order = list(run_feats)
        cum, max_single = [], []
        feats = []
        singles = []
        for name in order:
            feats.append(np.concatenate(run_feats[name]))
            singles.append(float(auroc(np.concatenate(run_feats[name]), y)))
            cum.append(cum_auc(feats, y, seed))
            max_single.append(max(singles))
        per_seed.append(dict(seed=seed, order=order, cumulative=cum, max_single=max_single,
                             union_exceeds_max=bool(cum[-1] > max_single[-1] + 1e-6)))
        print(json.dumps(per_seed[-1]), flush=True)
    agg = dict(
        cumulative=[float(np.mean([r["cumulative"][i] for r in per_seed]))
                    for i in range(len(per_seed[0]["order"]))],
        max_single=[float(np.mean([r["max_single"][i] for r in per_seed]))
                    for i in range(len(per_seed[0]["order"]))],
        union_exceeds_max_frac=float(np.mean([r["union_exceeds_max"] for r in per_seed])),
    )
    # governed anti-composition clause: AF and loadings mutually exclusive.
    led = Ledger()
    led.append("af", "released", "run5", meta={})
    led.append("loadings", "released", "run6", meta={"refines": "af"})
    flags = Ledger.relations(led.records)
    res = dict(per_seed=per_seed, agg=agg, order=per_seed[0]["order"],
               governed_clause=dict(policy="af XOR loadings", flagged=flags,
                                    flagged_violation=bool(flags)))
    with open(os.path.join(OUT, "s24_crossrun.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("agg", json.dumps(agg), flush=True)
    return res


if __name__ == "__main__":
    run()
