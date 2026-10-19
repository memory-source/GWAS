"""S2c - Mechanism-corrected marginal leakage of the PCA-diagnostic release.

The S2b delta comes from the released PCA being *study-fitted* - members'
projections onto their own fitted loadings sit systematically FARTHER from the
released centroid than non-members'.
Controls reported:
  T1  per-feature AUC (base / pcdist / pcalign / M2-only)
  T2  independent-spectrum control (loadings fit on attacker ref; both target
      panels out-of-sample) - isolates the overfitting component
  T3  released-spectrum attack (study-fit loadings + public study mean) - what
      an attacker can actually compute from the release
Output: results/raw/s2c_mechanism.json
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
    return {"pcdist": -np.linalg.norm(proj, axis=1), "pcalign": np.abs(proj[:, 0])}


def run(seeds=(0, 1, 2, 3, 4)):
    out = {"per_seed": []}
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, nonmem = split_population3(pop, 300, 300, 300)
        rel = release_objects(site)
        mask = rel["mask"]
        r_af = allele_freq(ref["G"])
        s_af = np.round(rel["af"], 2)
        y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(nonmem["G"]))])
        study_mean = site["G"][:, mask].mean(0)

        # The study-fitted PC basis is part of the released object (federated
        # PCA shares loadings so sites can project their own genotypes). All
        # study-fitted attacks below consume rel["loadings"], i.e. the released
        # object -- not a fresh recomputation from raw genotypes.
        eig_s, load_s = rel["eig"], rel["loadings"]
        _, _, load_r = pca(ref["G"], k=10)

        def base(P):
            return af_lrt_attack(P["G"], s_af, r_af, mask=mask)

        def pool(dm, dn, keys=None):
            if keys is None:
                keys = list(dm)
            return np.column_stack([np.concatenate([dm[k], dn[k]]) for k in keys])

        bm, bn = {"af": base(site)}, {"af": base(nonmem)}
        # T1: per-feature (s2b-style, per-panel centering)
        fmm = dict(bm); fmm.update(m2(site["G"], mask, load_s, site["G"][:, mask].mean(0)))
        fnn = dict(bn); fnn.update(m2(nonmem["G"], mask, load_s, nonmem["G"][:, mask].mean(0)))
        auc_base = cv_auc(pool(bm, bn, ["af"]), y, seed)
        auc_pcd = cv_auc(pool({"d": fmm["pcdist"]}, {"d": fnn["pcdist"]}), y, seed)
        auc_pca = cv_auc(pool({"a": fmm["pcalign"]}, {"a": fnn["pcalign"]}), y, seed)
        auc_m2only = cv_auc(pool({"d": fmm["pcdist"], "a": fmm["pcalign"]},
                                 {"d": fnn["pcdist"], "a": fnn["pcalign"]}), y, seed)
        auc_baseM2 = cv_auc(pool(fmm, fnn, ["af", "pcdist", "pcalign"]), y, seed)
        # T2: independent-spectrum control
        gmm = dict(bm); gmm.update(m2(site["G"], mask, load_r, study_mean))
        gnn = dict(bn); gnn.update(m2(nonmem["G"], mask, load_r, study_mean))
        auc_T2 = cv_auc(pool(gmm, gnn, ["af", "pcdist", "pcalign"]), y, seed)
        # T3: released-spectrum (what attacker computes)
        hmm = dict(bm); hmm.update(m2(site["G"], mask, load_s, study_mean))
        hnn = dict(bn); hnn.update(m2(nonmem["G"], mask, load_s, study_mean))
        auc_T3 = cv_auc(pool(hmm, hnn, ["af", "pcdist", "pcalign"]), y, seed)
        eig_ratio = float(eig_s.sum() / pca(ref["G"], k=10)[0].sum())
        # T5: within-ancestry label-permutation sanity (expect ~0.5)
        anc_all = np.concatenate([site["ancestry"], nonmem["ancestry"]])
        rng = np.random.default_rng(seed)
        perm = y.copy()
        for s in np.unique(anc_all):
            idx = np.where(anc_all == s)[0]
            perm[idx] = rng.permutation(perm[idx])
        auc_perm = cv_auc(pool(fmm, fnn, ["af", "pcdist", "pcalign"]), perm, seed)

        row = dict(seed=seed, base=auc_base, pcdist_only=auc_pcd, pcalign_only=auc_pca,
                   m2_only=auc_m2only, base_M2=auc_baseM2,
                   delta_M2=auc_baseM2 - auc_base,
                   independent_spectrum_T2=auc_T2, delta_T2=auc_T2 - auc_base,
                   released_spectrum_T3=auc_T3, delta_T3=auc_T3 - auc_base,
                   permuted_label_auc=auc_perm, eig_ratio=eig_ratio)
        out["per_seed"].append(row)
        print(json.dumps(row), flush=True)
    for k in ("base", "pcdist_only", "pcalign_only", "m2_only", "base_M2",
              "delta_M2", "independent_spectrum_T2", "delta_T2",
              "released_spectrum_T3", "delta_T3", "permuted_label_auc", "eig_ratio"):
        out[k + "_mean"] = float(np.mean([r[k] for r in out["per_seed"]]))
        out[k + "_std"] = float(np.std([r[k] for r in out["per_seed"]]))
    with open(os.path.join(OUT, "s2c_mechanism.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out


if __name__ == "__main__":
    run()
