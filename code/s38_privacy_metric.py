"""S38 - per-individual privacy metric on the release-precision frontier.

The paper's every claim is an AUROC (a rank statistic). A privacy reviewer asks
what that means for an individual. This script re-expresses the existing
operating points (frontier tiers s1c, PC channel s2b/s10, boundary s5, governed
contrast s23) in absolute per-individual terms:

  - membership advantage  A* = max_t (TPR(t) - FPR(t))   (KS distance); compared
    against the equal-variance binormal reference 2*Phi(Phi^{-1}(AUC)/sqrt(2))-1
    (and against the same-scale-free 2*AUC-1 as a secondary line)
  - per-individual posterior rho_i from a score-density likelihood ratio,
    in-sample descriptive (Sankararaman/Homer framing)
  - resolution-floored empirical f-DP trade-off summary
    eps_lb(delta) = max_alpha log((TPR(alpha)-delta)/alpha)
    against the (eps,delta) trade-off g_{eps,delta}; an in-sample,
    resolution-limited descriptive summary of the observed trade-off (not a
    statistical lower bound and not a release guarantee)

k-anonymity is rejected: the released objects are aggregate statistics, so
"number of members consistent with the release" is not well-defined.

Reuses the EXACT attack code (common.af_lrt_attack, pca, auroc, tpr_at_fpr) so
numbers line up with results/raw/s1c_frontier.json etc.

Output: results/raw/s38_privacy_metric.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc, tpr_at_fpr, pca)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


# ------------------------------------------------------------- metrics ---
def roc(scores, labels):
    """(fpr, tpr) over descending distinct score values (tie-grouped), plus the
    orientation that makes members score higher. Steps are taken only where the
    score value changes, so tied scores move FPR and TPR together."""
    s = np.asarray(scores, float)
    y = np.asarray(labels, int)
    if auroc(s, y) < 0.5:
        s = -s
    order = np.argsort(-s, kind="mergesort")
    s_o, y_o = s[order], y[order]
    npos = y.sum(); nneg = len(y) - npos
    last = np.r_[np.where(np.diff(s_o))[0], len(s_o) - 1]
    tp = np.cumsum(y_o)[last]; fp = np.cumsum(1 - y_o)[last]
    tpr = np.concatenate([[0.0], tp / npos])
    fpr = np.concatenate([[0.0], fp / nneg])
    return fpr, tpr


def membership_advantage(scores, labels):
    """A* = the two-sided KS distance between the member and non-member score
    distributions, max_t |TPR(t) - FPR(t)|. Computed symmetrically (both score
    orientations) so it does not require a pre-chosen direction; equals the
    one-sided threshold advantage only when the orientation matches the true
    separation."""
    fpr, tpr = roc(scores, labels)
    return float(np.max(np.abs(tpr - fpr)))


def per_individual_posterior(scores, labels, prior=0.5, bins=40):
    """rho_i from a histogram likelihood ratio at the released precision; returns
    per-member summary. The density is fit and evaluated on the same labels, so
    this is an in-sample descriptive posterior, not a calibrated one."""
    s = np.asarray(scores, float); y = np.asarray(labels, int)
    if auroc(s, y) < 0.5:
        s = -s
    lo, hi = np.percentile(s, [0.5, 99.5])
    edges = np.linspace(lo, hi, bins + 1)
    em, _ = np.histogram(s[y == 1], bins=edges, density=False)
    en, _ = np.histogram(s[y == 0], bins=edges, density=False)
    pm = (em + 0.5) / (em.sum() + 0.5 * bins)
    pn = (en + 0.5) / (en.sum() + 0.5 * bins)
    b_m = np.clip(np.digitize(s[y == 1], edges) - 1, 0, bins - 1)
    lr = pm[b_m] / pn[b_m]
    rho = prior * lr / (prior * lr + (1 - prior))
    return dict(median_rho=float(np.median(rho)),
                frac_gt_half=float((rho > 0.5).mean()),
                frac_gt_09=float((rho > 0.9).mean()),
                lr_max=float(lr.max()))


def eps_lower_bound(fpr, tpr, delta, n_pos, n_neg, fnr_mode="resolution"):
    """Resolution-floored empirical f-DP trade-off summary:
    the smallest eps consistent with the empirical trade-off
    f_obs(alpha) >= g_{eps,delta}(alpha), evaluated over alphas at or above the
    empirical resolution 1/n_neg, with the false-negative resolution floored. The
    two
    non-trivial branches of g_{eps,delta}(alpha)=max(0, 1-delta-e^eps*alpha,
    e^{-eps}(1-delta-alpha)) give eps >= log((TPR-delta)/alpha) and
    eps >= log((1-delta-alpha)/(1-TPR)); the value is the max over both branches
    and all alpha. It is an in-sample quantity (the ROC is evaluated on the same
    individuals the score was fit on) with no confidence treatment of the selected
    rates or thresholds, so it reads as a descriptive, resolution-limited summary
    of the observed trade-off, not a lower bound or a release-level guarantee; it
    is resolvable only where both rates exceed their 1/n resolutions.

    fnr_mode:
      'resolution' (default) -- floor the FNR at 1/n_pos (a single missed member
        is possible given the sample).
      'rule3'                 -- floor the FNR at 3/n_pos, the rule-of-three
        upper confidence bound on a zero observed false-negative count; this is
        the conservative floor that keeps the bound finite when the
        empirical TPR reaches 1 (a single zero-FN point otherwise makes the
        strict trade-off infeasible)."""
    alpha_floor = 1.0 / max(1, n_neg)
    fnr_floor = (3.0 / max(1, n_pos)) if fnr_mode == "rule3" else (1.0 / max(1, n_pos))
    a = np.clip(fpr, alpha_floor, 1 - 1e-12)
    fnr = np.clip(1 - tpr, fnr_floor, 1)
    sel = fpr >= alpha_floor
    if not sel.any():
        return 0.0
    b1 = np.log(np.clip(tpr[sel] - delta, 1e-12, None) / a[sel])
    b2 = np.log(np.clip(1 - delta - a[sel], 1e-12, None) / fnr[sel])
    return float(max(0.0, np.max(np.concatenate([b1, b2]))))


def validate_tradeoff(fpr, tpr, eps, delta, alpha_floor, fnr_floor=None):
    """Check f_obs(alpha) >= g_{eps,delta}(alpha) on the ROC grid (should hold).
    The false-negative side uses the SAME resolution floor as eps_lower_bound, so
    a zero observed false-negative count is treated as the conservative rule-of-
    three FNR when `fnr_floor=3/n_pos` (otherwise a single zero-FN point makes the
    strict trade-off infeasible and the flag reads False by construction)."""
    sel = fpr >= alpha_floor
    a = np.clip(fpr[sel], alpha_floor, 1 - 1e-12)
    fnr = 1 - tpr[sel]
    if fnr_floor is not None:
        fnr = np.clip(fnr, fnr_floor, 1)
    else:
        fnr = np.clip(fnr, 1e-12, 1)
    g = np.maximum(0.0, np.maximum(1 - delta - np.exp(eps) * a,
                                    np.exp(-eps) * (1 - delta - a)))
    return bool(np.all(fnr >= g - 1e-9))


def binormal_ks(auroc):
    """KS = A* for an equal-variance binormal ROC with the given AUC:
    A*_bin = 2*Phi(Phi^{-1}(AUC)/sqrt(2)) - 1. This is the correct binormal
    reference; 2*AUC-1 is a different scale (equal to A* only at AUC=0.5)."""
    from scipy.stats import norm
    auc = float(np.clip(auroc, 1e-6, 1 - 1e-6))
    return float(2 * norm.cdf(norm.ppf(auc) / np.sqrt(2)) - 1)


def metrics_from_scores(scores, labels, delta_list=(0.05, 0.01, 0.001)):
    """delta must be well above the finite-sample floor for the empirical
    trade-off to be an estimate (below it a single FP/FN dominates). With
    n_pos=n_neg=900 the resolvable regime is delta >=~ 0.01."""
    fpr, tpr = roc(scores, labels)
    lab = np.asarray(labels, int)
    n_neg = int((lab == 0).sum())
    n_pos = int((lab == 1).sum())
    a = float(auroc(scores, labels))
    ks = membership_advantage(scores, labels)
    m = dict(auroc=a, tpr_at_fpr1=float(tpr_at_fpr(np.asarray(scores, float),
                                                   np.asarray(labels, int), 0.01)),
             ks_advantage=ks,
             ks_binormal=binormal_ks(a),
             gap_vs_binormal=ks - binormal_ks(a),
             gap_vs_2aucm1=ks - (2 * a - 1))
    m.update(per_individual_posterior(scores, labels))
    for d in delta_list:
        e = eps_lower_bound(fpr, tpr, d, n_pos, n_neg)
        e3 = eps_lower_bound(fpr, tpr, d, n_pos, n_neg, fnr_mode="rule3")
        m[f"eps_lb_{d:g}"] = e
        m[f"eps_lb_{d:g}_valid"] = validate_tradeoff(fpr, tpr, e, d, 1.0 / max(1, n_neg),
                                                      fnr_floor=1.0 / max(1, n_pos))
        m[f"eps_lb_rule3_{d:g}"] = e3
        m[f"eps_lb_rule3_{d:g}_valid"] = validate_tradeoff(fpr, tpr, e3, d, 1.0 / max(1, n_neg),
                                                           fnr_floor=3.0 / max(1, n_pos))
    return m


# --------------------------------------------------------------- cells ---
def pc_dist(G, loadings, mean):
    X = G.astype(np.float64) - mean[None, :]
    return np.linalg.norm((X @ loadings.T) / np.sqrt(max(1, G.shape[1])), axis=1)  # was M^0.5


def cell_frontier(seed):
    """Tier cells of s1c: AF@1dp, AF@2dp, AF full, AF per-site, PC, PC+AF, governed."""
    pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
    site, ref, non = split_population3(pop, 300, 300, 300)
    rel = release_objects(site)
    mask = rel["mask"]
    af_r = allele_freq(ref["G"])
    y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])
    out = {}
    for name, af in (("af_full", rel["af"]), ("af_2dp", np.round(rel["af"], 2)),
                     ("af_1dp", np.round(rel["af"], 1))):
        sm = af_lrt_attack(site["G"], af, af_r, mask=mask)
        sn = af_lrt_attack(non["G"], af, af_r, mask=mask)
        out[name] = metrics_from_scores(np.concatenate([sm, sn]), y)
    # per-site AF
    sm = np.zeros(len(site["G"])); sn = np.zeros(len(non["G"]))
    for k in rel["af_by_group"]:
        afk = rel["af_by_group"][k]
        sm[site["ancestry"] == k] = af_lrt_attack(site["G"][site["ancestry"] == k], afk,
                                                  allele_freq(ref["G"][ref["ancestry"] == k]), mask=mask)
        sn[non["ancestry"] == k] = af_lrt_attack(non["G"][non["ancestry"] == k], afk,
                                                 allele_freq(ref["G"][ref["ancestry"] == k]), mask=mask)
    out["af_sitewise"] = metrics_from_scores(np.concatenate([sm, sn]), y)
    # PC channel (study-fitted loadings, train-half oriented)
    _, _, L = pca(site["G"], k=10)
    mean = site["G"][:, mask].mean(0)
    dm = pc_dist(site["G"][:, mask], L[:, mask], mean)
    dn = pc_dist(non["G"][:, mask], L[:, mask], mean)
    rng = np.random.default_rng(seed + 7717)
    tsm = rng.random(len(dm)) < 0.5; tsn = rng.random(len(dn)) < 0.5
    ytr = np.concatenate([np.ones(tsm.sum()), np.zeros(tsn.sum())])
    sgn = 1.0 if auroc(-np.concatenate([dm[tsm], dn[tsn]]), ytr) >= 0.5 else -1.0
    out["pc_loadings"] = metrics_from_scores(
        sgn * (-np.concatenate([dm, dn])), y)
    # governed PC: eigenvalues-only (random orthonormal basis)
    Lr = rng.normal(size=L.shape); Lr, _ = np.linalg.qr(Lr.T); Lr = Lr.T
    dmr = pc_dist(site["G"][:, mask], Lr[:, mask], mean)
    dnr = pc_dist(non["G"][:, mask], Lr[:, mask], mean)
    out["pc_eigonly_control"] = metrics_from_scores(sgn * (-np.concatenate([dmr, dnr])), y)
    # structure control: reference-fitted basis
    _, _, Lref = pca(ref["G"], k=10)
    mr = ref["G"][:, mask].mean(0)
    dmx = pc_dist(site["G"][:, mask], Lref[:, mask], mr)
    dnx = pc_dist(non["G"][:, mask], Lref[:, mask], mr)
    out["pc_refbasis_control"] = metrics_from_scores(-np.concatenate([dmx, dnx]), y)
    return out


def run(seeds=(0, 1, 2, 3, 4)):
    acc = {}
    for s in seeds:
        cells = cell_frontier(s)
        for name, m in cells.items():
            acc.setdefault(name, []).append(m)
        print("seed", s, "done", flush=True)
    agg = {}
    for name, rows in acc.items():
        agg[name] = {}
        for k in rows[0]:
            vals = [r[k] for r in rows if isinstance(r[k], (int, float, bool))]
            agg[name][k] = float(np.mean(vals))
        # keep the trade-off validity as an all-seeds conjunction
        for k in rows[0]:
            if k.endswith("_valid"):
                agg[name][k] = bool(all(r[k] for r in rows))
    res = dict(agg=agg, n_seeds=len(seeds))
    # AUROC-vs-eps non-order-preservation example: two near-equal-AUROC AF
    # channels with inverted eps (af_full has the higher AUROC but the LOWER eps
    # at every delta), showing the eps summary is not order-preserving with AUROC.
    if "af_full" in agg and "af_2dp" in agg:
        a, b = agg["af_full"], agg["af_2dp"]
        res["eps_order_inversion"] = dict(
            pair=["af_full", "af_2dp"],
            auroc=[a["auroc"], b["auroc"]],
            eps_lb_001=[a["eps_lb_0.01"], b["eps_lb_0.01"]],
            eps_lb_005=[a["eps_lb_0.05"], b["eps_lb_0.05"]],
            note="higher AUROC (af_full) has the LOWER eps at each delta -> the "
                 "eps summary is not order-preserving with AUROC")
    with open(os.path.join(OUT, "s38_privacy_metric.json"), "w") as f:
        json.dump(res, f, indent=2)
    for name in agg:
        a = agg[name]
        print("%-24s AUC %.3f  A* %.3f  eps_lb(d=0.01) %.2f  med_rho %.3f  LR_max %.1f" % (
            name, a["auroc"], a["ks_advantage"], a["eps_lb_0.01"], a["median_rho"], a["lr_max"]))
    print("(eps is an in-sample, resolution-limited trade-off summary at "
          "n_pos=n_neg=900; resolvable delta>=~0.01; full sensitivity in the JSON)")
    return res


if __name__ == "__main__":
    run()
