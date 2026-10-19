"""S10 - paired significance for the PC channel vs the AF-only channel.

Result-to-claim gap #2: the PC-vs-AF comparison was reported as an
unpaired mean gap (+0.037 across 5 seeds) with no formal test. Here we run
N>=20 seeds, pair every channel within a seed, calibrate the PC orientation
on a held-out TRAIN half (so the reported test AUROC uses no test labels -
this removes the "oracle orientation" objection), and report paired tests
(t-test, Wilcoxon, sign test) plus a bootstrap CI on the gap.

Channels (per sample):
  af_only       : Homer/Sankararaman AF-LRT, reference-panel frequencies
  pc_calibrated : study-fitted projection distance (full-study basis), sign
                  calibrated on the train half
  combined      : logistic regression on [af, pc_calibrated] fit on train
  af_2dp        : AF-LRT using 2-dp-rounded released AF
  af_conditioned: AF-LRT against the study AF with a composition-balanced subset
                  of the KNOWN train members' genotypes subtracted; the residual
                  frequency is used for BOTH classes (disjoint negatives). A
                  valid conditioned-AF attack.
  pc_conditioned: the study PC basis REFIT on the SAME remaining panel and the
                  same projection-distance attack. This is a smaller-cohort
                  CURATOR-refitted release experiment (it changes the released
                  object using the remaining members' genotypes, which are not
                  recoverable from the original released statistics), NOT an
                  attacker reconstruction of the original PCA release. It is
                  reported separately from the as-released comparison; the
                  orientation sign is reused from the training-half calibration.
  af_loo_oracle : member-only leave-own-out residual (the scored member's own
                  genotype subtracted before scoring). ORACLE ablation: choosing
                  that branch uses the membership being predicted, so it is a
                  mechanism reading, not a deployable attack.

Output: results/raw/s10_paired.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc, pca, tpr_at_fpr)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def lr_fit_score(Xtr, ytr, Xte):
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler().fit(Xtr)
    clf = LogisticRegression(max_iter=2000).fit(sc.transform(Xtr), ytr)
    return clf.predict_proba(sc.transform(Xte))[:, 1]


def one_seed(seed):
    pop = simulate_population(n_groups=3, n_per_group=1200, n_snps=4000, seed=seed)
    site, ref, non = split_population3(pop, 300, 300, 300)
    rel = release_objects(site)
    mask = rel["mask"]
    af_s, af_r = rel["af"], allele_freq(ref["G"])
    mean_s = site["G"][:, mask].mean(0)
    load_s = rel["loadings"]

    def pc_score(G):
        X = G[:, mask].astype(np.float64) - mean_s[None, :]
        proj = (X @ load_s[:, mask].T) / np.sqrt(max(1, mask.sum()))
        return np.linalg.norm(proj, axis=1)

    dm, dn = pc_score(site["G"]), pc_score(non["G"])
    am = af_lrt_attack(site["G"], af_s, af_r, mask=mask)
    an = af_lrt_attack(non["G"], af_s, af_r, mask=mask)
    bm2 = af_lrt_attack(site["G"], np.round(af_s, 2), af_r, mask=mask)
    bn2 = af_lrt_attack(non["G"], np.round(af_s, 2), af_r, mask=mask)

    # deterministic train/test split of membership (balanced within each panel)
    rng = np.random.default_rng(seed + 999)
    tsm = rng.random(len(dm)) < 0.5
    tsn = rng.random(len(dn)) < 0.5
    ytr = np.concatenate([np.ones(tsm.sum()), np.zeros(tsn.sum())])
    yte = np.concatenate([np.ones((~tsm).sum()), np.zeros((~tsn).sum())])

    af_tr = np.concatenate([am[tsm], an[tsn]])
    af_te = np.concatenate([am[~tsm], an[~tsn]])
    pc_tr = np.concatenate([dm[tsm], dn[tsn]])
    pc_te = np.concatenate([dm[~tsm], dn[~tsn]])
    pc2_te = np.concatenate([bm2[~tsm], bn2[~tsn]])

    # calibrated-PC orientation on train only
    a_naive = auroc(pc_tr, ytr)
    a_flip = auroc(-pc_tr, ytr)
    sgn = 1.0 if a_naive >= a_flip else -1.0
    pc_cal_tr = sgn * pc_tr
    pc_cal_te = sgn * pc_te

    # Conditioned channels. The calibrated PC/logistic channels use the known
    # membership of the train half, so we give the AF channel the same
    # known-member side information and also report a curator-refitted PC arm:
    #   af_conditioned : subtract a composition-balanced subset of the known
    #                    train members from the RELEASED study AF and score the
    #                    residual frequency for both classes (one frequency,
    #                    disjoint negatives). This is attacker post-processing of
    #                    the original release.
    #   pc_conditioned : the study PC basis REFIT on the remaining members, then
    #                    the same projection-distance attack. This is a
    #                    smaller-cohort CURATOR-refitted release experiment (it
    #                    changes the released object using the remaining members'
    #                    genotypes, which are not recoverable from the original
    #                    released statistics), NOT an attacker reconstruction of
    #                    the original PCA release.
    anc = site["ancestry"]
    N = site["G"].shape[0]
    cnt_tr = [int((tsm & (anc == g)).sum()) for g in np.unique(anc)]
    bal = np.concatenate([np.where((tsm) & (anc == g))[0][:min(cnt_tr)]
                          for g in np.unique(anc)])
    rem_mask = np.zeros(N, bool); rem_mask[bal] = True
    num2 = 2.0 * N * af_s - site["G"][bal].sum(0).astype(float)
    af_cond = np.clip(num2 / (2.0 * (N - len(bal))), 1e-6, 1 - 1e-6)
    cm_cond = af_lrt_attack(site["G"], af_cond, af_r, mask=mask)
    cn_cond = af_lrt_attack(non["G"], af_cond, af_r, mask=mask)
    cond_te = np.concatenate([cm_cond[~tsm], cn_cond[~tsn]])

    # PC refit on the remaining (kept) members -- a curator-refitted
    # smaller-cohort release (not an attack on the original release)
    kept = ~rem_mask
    mean_k = site["G"][kept][:, mask].mean(0)
    _, _, load_k = pca(site["G"][kept], k=10)

    def pc_cond_score(Gp):
        X = Gp[:, mask].astype(np.float64) - mean_k[None, :]
        proj = (X @ load_k[:, mask].T) / np.sqrt(max(1, mask.sum()))
        return np.linalg.norm(proj, axis=1)
    dkm, dkn = pc_cond_score(site["G"]), pc_cond_score(non["G"])
    # orientation: the balanced removal consumes almost all train members, so a
    # train-half re-orientation on the kept panel is noise (10 members). Reuse
    # the training-calibrated sign `sgn` fixed by the full-study train half.
    pc_cond_te = sgn * np.concatenate([dkm[~tsm], dkn[~tsn]])

    comb_te = lr_fit_score(np.column_stack([af_tr, pc_cal_tr]), ytr,
                           np.column_stack([af_te, pc_cal_te]))

    # oracle leave-own-out residual: recompute af_cond with the scored member's
    # own genotype ALSO subtracted. This uses the evaluation membership to choose
    # the branch (members are scored leave-own-out, non-members are not), so it is
    # a mechanism ablation, not a deployable attack; retained only to show what
    # self-contribution contributes.
    gK_all = site["G"].astype(np.float64)
    m_idx = np.where(mask)[0]
    loo_scores = np.empty(N)
    kept_n = N - len(bal)
    for i in range(N):
        num_i = num2[m_idx] - gK_all[i, m_idx]
        af_i = np.clip(num_i / (2.0 * (kept_n - 1)), 1e-6, 1 - 1e-6)
        g = gK_all[i, m_idx]
        lp = np.where(g == 0, np.log((1 - af_i) ** 2),
                      np.where(g == 1, np.log(2 * af_i * (1 - af_i)),
                               np.log(af_i ** 2)))
        lr = np.where(g == 0, np.log((1 - af_r[m_idx]) ** 2),
                      np.where(g == 1, np.log(2 * af_r[m_idx] * (1 - af_r[m_idx])),
                               np.log(af_r[m_idx] ** 2)))
        loo_scores[i] = (lp - lr).sum()
    af_loo_oracle = float(auroc(np.concatenate([loo_scores[~tsm], cn_cond[~tsn]]), yte))

    return dict(
        seed=seed,
        af_only=float(auroc(af_te, yte)),
        af_2dp=float(auroc(pc2_te, yte)),
        pc_calibrated=float(auroc(pc_cal_te, yte)),
        af_conditioned=float(auroc(cond_te, yte)),
        pc_conditioned=float(auroc(pc_cond_te, yte)),
        af_loo_oracle=af_loo_oracle,
        combined=float(auroc(comb_te, yte)),
        pc_sign=float(sgn),
        train_pc_naive=float(a_naive),
        train_pc_flip=float(a_flip),
        n_train_members=int(len(bal)),
    )


def paired_tests(a, b, seed=0):
    """Return paired-difference summaries + t-test / Wilcoxon / sign-test."""
    from scipy import stats
    a = np.asarray(a, float); b = np.asarray(b, float)
    d = a - b
    n = len(d)
    out = dict(mean=float(d.mean()), std=float(d.std(ddof=1)),
               n=n, mean_a=float(a.mean()), mean_b=float(b.mean()))
    out["t_p"] = float(stats.ttest_rel(a, b).pvalue)
    try:
        out["wilcoxon_p"] = float(stats.wilcoxon(a, b).pvalue)
    except Exception:
        out["wilcoxon_p"] = float("nan")
    out["sign_test_p"] = float(stats.binomtest(int((d > 0).sum()), n, 0.5).pvalue)
    rng = np.random.default_rng(seed)
    boot = [d[rng.integers(0, n, n)].mean() for _ in range(5000)]
    out["gap_ci95"] = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]
    return out


def run(n_seeds=25):
    rows = [one_seed(s) for s in range(n_seeds)]
    for r in rows:
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v)
                          for k, v in r.items()}), flush=True)
    chans = {c: [r[c] for r in rows] for c in
             ("af_only", "af_2dp", "pc_calibrated", "af_conditioned",
              "pc_conditioned", "af_loo_oracle", "combined")}
    res = dict(per_seed=rows, n_seeds=n_seeds,
               means={c: float(np.mean(v)) for c, v in chans.items()},
               stds={c: float(np.std(v)) for c, v in chans.items()})
    res["pc_vs_af"] = paired_tests(chans["pc_calibrated"], chans["af_only"])
    res["pc_vs_af2dp"] = paired_tests(chans["pc_calibrated"], chans["af_2dp"])
    res["pc_vs_afcond"] = paired_tests(chans["pc_calibrated"], chans["af_conditioned"])
    res["pccond_vs_afcond"] = paired_tests(chans["pc_conditioned"], chans["af_conditioned"])
    res["combined_vs_af"] = paired_tests(chans["combined"], chans["af_only"])
    res["pc_beats_afcond_frac"] = float(np.mean(
        [r["pc_calibrated"] > r["af_conditioned"] for r in rows]))
    res["pccond_beats_afcond_frac"] = float(np.mean(
        [r["pc_conditioned"] > r["af_conditioned"] for r in rows]))
    res["sign_agreement"] = float(np.mean([r["pc_sign"] > 0 for r in rows]))
    res["mean_train_pc_naive"] = float(np.mean([r["train_pc_naive"] for r in rows]))
    with open(os.path.join(OUT, "s10_paired.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("=== summary ===", flush=True)
    for k in ("means", "stds", "pc_vs_af", "combined_vs_af", "sign_agreement"):
        print(k, json.dumps(res[k] if not isinstance(res[k], dict)
                            else {kk: (round(vv, 4) if isinstance(vv, float) else vv)
                                  for kk, vv in res[k].items()}), flush=True)
    return res


if __name__ == "__main__":
    run()
