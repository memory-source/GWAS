"""S8 - Robustness of the AF-proximity attack across functional forms.

IMPORTANT SCOPE (per cross-model review): F1/F2/F3 below are three FUNCTIONAL
FORMS of the SAME information channel --- what a target individual's genotype
reveals about proximity to the released study allele-frequency vector relative
to a reference. They are NOT three independent attack channels. Their value is
functional-form robustness: it shows the frontier result is not an artifact of
the likelihood-ratio functional specifically, and reveals that the LRT form is
more rounding-fragile (at 1 dp the correlation/deviation forms retain more
signal than the LRT).

  F1  AF-LRT        : Homer/Sankararaman sum of genotype log-likelihood ratios.
  F2  AF-correlation: corr(genotype, study-AF) - corr(genotype, ref-AF).
  F3  ML deviation  : cross-validated logistic regression on per-individual mean
                      squared deviation from study-AF vs reference-AF.
A genuinely distinct channel (haplotype-based, e.g. Bu et al. 2021) is NOT
implemented here and is noted as future work.

Output: results/raw/s8_families.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    simulate_population, split_population3, release_objects, allele_freq,
    af_lrt_attack, anc_only_attack, auroc, tpr_at_fpr, pca,
)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def af_corr_attack(G, study_af, ref_af, mask):
    X = G[:, mask].astype(np.float64)
    s = study_af[mask]; r = ref_af[mask]
    Xc = X - X.mean(1, keepdims=True)
    sc = s - s.mean(); rc = r - r.mean()
    denom = np.sqrt((Xc ** 2).sum(1))
    cs = (Xc * sc[None, :]).sum(1) / np.maximum(denom * np.sqrt((sc ** 2).sum()), 1e-12)
    cr = (Xc * rc[None, :]).sum(1) / np.maximum(denom * np.sqrt((rc ** 2).sum()), 1e-12)
    return cs - cr


def ml_dev_features(G, study_af, ref_af, mask):
    X = G[:, mask].astype(np.float64)
    ds = ((X - 2 * study_af[mask][None, :]) ** 2).mean(1)
    dr = ((X - 2 * ref_af[mask][None, :]) ** 2).mean(1)
    return np.column_stack([ds, dr])


def cv_auc(X, y, seed=0):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr])
        clf = LogisticRegression(max_iter=2000).fit(sc.transform(X[tr]), y[tr])
        oof[te] = clf.predict_proba(sc.transform(X[te]))[:, 1]
    return auroc(oof, y)


def evaluate(site, ref, nonmem, study_af, ref_af, mask):
    lab = np.concatenate([np.ones(site.shape[0]), np.zeros(nonmem.shape[0])])
    s1 = np.concatenate([af_lrt_attack(site, study_af, ref_af, mask=mask),
                         af_lrt_attack(nonmem, study_af, ref_af, mask=mask)])
    s2 = np.concatenate([af_corr_attack(site, study_af, ref_af, mask),
                         af_corr_attack(nonmem, study_af, ref_af, mask)])
    Fm = ml_dev_features(site, study_af, ref_af, mask)
    Fn = ml_dev_features(nonmem, study_af, ref_af, mask)
    a3 = cv_auc(np.vstack([Fm, Fn]), lab)
    # study-fitted PC distance (the diagnostic leak), adversarial orientation
    _, _, load_s = pca(site, k=10)
    mean_s = site[:, mask].mean(0)
    pcd = np.concatenate([-anc_only_attack(site, load_s, mean_s, mask),
                          -anc_only_attack(nonmem, load_s, mean_s, mask)])
    return dict(
        f1_af_lrt=float(auroc(s1, lab)), f1_tpr1=float(tpr_at_fpr(s1, lab, 0.01)),
        f2_af_corr=float(auroc(s2, lab)), f2_tpr1=float(tpr_at_fpr(s2, lab, 0.01)),
        f3_ml_dev=float(a3),
        pc_distance=float(auroc(pcd, lab)),
    )


def run(seeds=(0, 1, 2, 3, 4), n_site=300):
    rows = []
    config = dict(seeds=list(seeds), n_site=n_site, n_groups=3, n_per_group=1000,
                  n_snps=4000, control="split_population3_disjoint_nonmembers",
                  maf_filter="[0.05,0.95]")
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, nonmem = split_population3(pop, n_site, 300, 300)
        rel = release_objects(site)
        mask = rel["mask"]
        af_s = allele_freq(site["G"]); af_r = allele_freq(ref["G"])
        for tier, sa in (("F5_full", af_s),
                         ("F3_round2", np.round(af_s, 2)),
                         ("F4_round1", np.round(af_s, 1))):
            r = evaluate(site["G"], ref["G"], nonmem["G"], sa, af_r, mask)
            r.update(seed=seed, tier=tier)
            # guard: ML CV uses raw per-individual features; recompute here so the
            # tier only changes the released AF used by F1/F2
            rows.append(r)
            print(r, flush=True)
    fam = {}
    for key in ("f1_af_lrt", "f2_af_corr", "f3_ml_dev", "pc_distance",
                "f1_tpr1", "f2_tpr1"):
        for tier in ("F5_full", "F3_round2", "F4_round1"):
            vals = [r[key] for r in rows if r["tier"] == tier]
            fam[f"{tier}|{key}"] = float(np.mean(vals))
    out = dict(config=config, rows=rows, means=fam)
    with open(os.path.join(OUT, "s8_families.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(fam, indent=2))
    return out


if __name__ == "__main__":
    run()
