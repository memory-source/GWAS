"""S67 - fair comparison against the FUSED passive adversary (round 22).

The round-19/20/21 margin was measured against the stronger of two SINGLE passive
channels (`max(af_cond, spec_cond)`). A real passive adversary holds both channels
and combines them, so the decisive comparison is the accumulated active attacker
against a FUSED passive adversary.

For each cohort x split seed:
  - build the two conditioned passive channels on the SAME held-out split as the
    spectral attack: `af_cond` (composition-balanced removal of the known
    calibration members from the released AF) and `spec_cond` (the conditioned
    single release with the calibration members' scatter removed);
  - sign each channel on the calibration half (its own unconditioned orientation);
  - FU SE: z-score both signed channels using the pooled row statistics (no labels),
    fused_z = z_af + z_spec; also fit a logistic fusion on the calibration half's
    labels (the attacker knows the calibration members);
  - the JOINT passive adversary is the candidate (each single channel, fused_z,
    logistic fusion) with the best CALIBRATION-half AUROC, evaluated held-out;
  - run the q = 1/2/4/8 random/adaptive stacks exactly as in s61.cell and score the
    membership from the stacked operator;
  - report `gap_vs_joint_passive` in AUC and TPR@5%FPR, with bootstrap CIs and the
    positive split count, plus the center/scale-only control (af_pooled) and the
    same-input passive-repeat control.

Output: results/raw/s67_joint_passive_chr<C>.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
import s51_real_lowstruct as s51
import s61_spectral_tomography as s61
from common import allele_freq, af_lrt_attack, auroc, tpr_at_fpr
from real_multi import SUPERPOP

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")


def cohorts_for(chrom, ld_window=1500, max_snps=20000):
    Ge, _, _ = s48.load_for(("EUR",), chrom, ld_window, max_snps)
    G3, _, sup3, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    G5, _, sup5, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    return {"within:EUR": (Ge, np.array(["EUR"] * len(Ge))),
            "super:3pop": (G3, sup3),
            "within:EAS": (G5[sup5 == "EAS"], sup5[sup5 == "EAS"])}


def _conditioned_af_raw(honest, nonmem, ref, partner, shadow, honest_groups, y, cal,
                        seed):
    """Raw conditioned-AF-LRT score for all (honest + nonmember) rows, plus the
    pooled-AF (center/scale-only) raw score and the unconditioned sign. Mirrors
    s61.af_baselines but returns raw values so they can be fused."""
    af_s = allele_freq(honest)
    bg_pooled = allele_freq(np.vstack([ref, partner, shadow]))
    h = np.asarray(cal, bool)
    raw_ref = np.concatenate([af_lrt_attack(honest, af_s, allele_freq(ref)),
                              af_lrt_attack(nonmem, af_s, allele_freq(ref))])
    raw_pooled = np.concatenate([af_lrt_attack(honest, af_s, bg_pooled),
                                 af_lrt_attack(nonmem, af_s, bg_pooled)])
    sgn_af = 1.0 if auroc(raw_pooled[h], y[h]) >= 0.5 else -1.0
    N = len(honest)
    gh = np.asarray(honest_groups)
    h_study = h[:N]
    groups_present = np.unique(gh)
    quota = min(int((h_study & (gh == g)).sum()) for g in groups_present)
    bal = np.concatenate([np.where(h_study & (gh == g))[0][:quota] for g in groups_present])
    num2 = 2.0 * N * af_s - honest[bal].sum(0).astype(float)
    af_cond = np.clip(num2 / (2.0 * (N - len(bal))), 1e-6, 1 - 1e-6)
    raw_cond = np.concatenate([af_lrt_attack(honest, af_cond, bg_pooled),
                               af_lrt_attack(nonmem, af_cond, bg_pooled)])
    return raw_cond, raw_pooled, sgn_af, int(len(bal))


def _projection_score(R, honest, nonmem, mu, scale):
    """Raw held-out score = projection norm onto span(R). R is M x d with d small
    (the released/stacked operator); the span is its column space, so the SVD is
    cheap. `sign` handled by the caller."""
    U, s, _ = np.linalg.svd(R, full_matrices=False)
    keep = s > 1e-9 * (s[0] if s.size else 1.0)
    B = U[:, keep]
    X = s61.zscore(np.vstack([honest, nonmem]), mu, scale)
    return np.linalg.norm(X @ B, axis=1)


def _z(x):
    sd = x.std()
    return (x - x.mean()) / (sd if sd > 0 else 1.0)


def one_split(G, groups, seed, k=s61.K_TOP, block_len=s61.BLOCK_LEN, caps=None):
    roles = s61.split_roles(groups, caps or s61.CAPS, seed)
    honest, partner = G[roles["honest"]], G[roles["partner"]]
    nonmem, ref, shadow = G[roles["nonmember"]], G[roles["reference"]], G[roles["shadow"]]
    mu, scale = s61.joint_baseline(honest, partner)
    Zh = s61.zscore(honest, mu, scale)
    Zsh = s61.zscore(shadow, mu, scale)
    n_blocks = G.shape[1] // block_len
    y = np.concatenate([np.ones(len(honest)), np.zeros(len(nonmem))])
    cal = np.random.default_rng(seed + 31).random(len(y)) < 0.5
    h = np.asarray(cal, bool)
    honest_groups = np.asarray(groups)[roles["honest"]]
    grp_rows = np.concatenate([honest_groups, np.asarray(groups)[roles["nonmember"]]])

    raw_cond, raw_pooled, sgn_af, n_bal = _conditioned_af_raw(
        honest, nonmem, ref, partner, shadow, honest_groups, y, cal, seed)
    s_af_cond = sgn_af * raw_cond
    s_pooled = sgn_af * raw_pooled

    # single passive release + its conditioned form
    Zp0 = s61.zscore(partner, mu, scale)
    lam0, V0 = s61.topk([Zh, Zp0], k)
    K0m = Zp0.T @ Zp0
    Y0 = V0 * lam0[None, :] - K0m @ V0
    Zk = Zh[h[:len(honest)]]
    Y0_cond = Y0 - Zk.T @ (Zk @ V0)
    _, _, sgn_spec = s61._membership_auc_from_operator(Y0, honest, nonmem, mu, scale, y, cal)
    s_spec = sgn_spec * _projection_score(Y0_cond, honest, nonmem, mu, scale)

    # ---- passive adversary candidates ----
    # The two CONDITIONED channels (af_cond, spec_cond) are the regime the reviewer
    # asked for. Each is oriented on its UNCONDITIONED calibration half (sgn_af /
    # sgn_spec, both legitimate; the conditioning removes the calibration members, so
    # the conditioned calibration half is uninformative). The adversary JOINS them.
    #   joint_passive  = fused_z, the label-free z-sum of the two conditioned channels
    #                    (weights are not fitted, so no calibration signal is needed)
    #   fused_cal      = a diagnostic fusion whose weights ARE fitted on the
    #                    unconditioned calibration half (af_pooled, spec_uncond) and
    #                    applied to the conditioned channels -- reported beside, not
    #                    as the primary baseline
    #   joint_oracle   = the best held-out single/joint candidate, an explicit oracle
    #                    upper bound (not reachable without held-out labels)
    def heldout(s):        # AUROC and TPR@5% on the held-out half
        return float(auroc(s[~h], y[~h])), float(tpr_at_fpr(s[~h], y[~h], 0.05))

    fused_z = _z(s_af_cond) + _z(s_spec)
    from sklearn.linear_model import LogisticRegression
    # diagnostic fusion: weights fitted on the STANDARDIZED unconditioned calibration
    # features and applied to the STANDARDIZED conditioned features (same units)
    s_spec_uncond = sgn_spec * _projection_score(Y0, honest, nonmem, mu, scale)
    Z_uncond = np.column_stack([_z(s_pooled), _z(s_spec_uncond)])
    w = LogisticRegression(max_iter=1000).fit(Z_uncond[h], y[h]).coef_[0]
    fused_cal = w[0] * _z(s_af_cond) + w[1] * _z(s_spec)

    joint_cands = {"af_cond": s_af_cond, "spec_cond": s_spec,
                   "fused_z": fused_z, "fused_cal": fused_cal}
    ho = {k: heldout(v) for k, v in joint_cands.items()}
    # PRIMARY AUC bar = the AUROC oracle: the candidate with the best HELD-OUT AUROC among
    # {af_cond, spec_cond, fused_z, fused_cal}. It is the strongest of these four by
    # AUROC, so no AUROC-selection rule can beat it (it is labeled an oracle because the
    # choice uses held-out labels). Its TPR is that candidate's TPR, not the max TPR of
    # the four. A metric-specific TPR oracle is reported separately (tpr_oracle_*), since
    # the comparator for the TPR reading should also be the strongest TPR passive bar.
    oracle_key = max(ho, key=lambda k: ho[k][0])
    joint_auc, joint_tpr = ho[oracle_key]
    tpr_oracle_key = max(ho, key=lambda k: ho[k][1])
    tpr_oracle = ho[tpr_oracle_key][1]
    tpr_fused_z = ho["fused_z"][1]
    tpr_fused_cal = ho["fused_cal"][1]

    # ---- active attacker: the q-stacks, scored from the stacked operator ----
    x_in, x_out, _, _ = s61._target_dirs(Zsh, np.linalg.qr(V0)[0])
    out = {"seed": seed, "n_conditioned_members": n_bal,
           "joint_candidates_auc": {k: ho[k][0] for k in ho},
           "joint_candidates_tpr05": {k: ho[k][1] for k in ho},
           "oracle_key": oracle_key,
           "tpr_oracle_key": tpr_oracle_key,
           "joint_passive_auc": joint_auc, "joint_passive_tpr05": joint_tpr,
           "tpr_oracle_bar": tpr_oracle,
           "fused_z_auc": ho["fused_z"][0], "fused_z_tpr05": tpr_fused_z,
           "fused_cal_auc": ho["fused_cal"][0], "fused_cal_tpr05": tpr_fused_cal,
           "af_cond_auc": heldout(s_af_cond)[0], "af_cond_tpr05": heldout(s_af_cond)[1],
           "spec_cond_auc": heldout(s_spec)[0], "spec_cond_tpr05": heldout(s_spec)[1],
           "pooled_center_scale_auc": heldout(s_pooled)[0],
           "pooled_center_scale_tpr05": heldout(s_pooled)[1],
           "strategies": {}}
    for strategy in ("passive", "random", "adaptive"):
        rng = np.random.default_rng(seed + 101)
        Qs, Ys = [V0], [Y0]
        per_q = {}
        for t in range(0, max(s61.QS)):
            if t > 0:
                if strategy == "passive":
                    own = partner
                elif strategy == "random":
                    own, _ = s61.swap_block(partner, rng.choice(partner.shape[0], 2,
                                                                replace=False),
                                            int(rng.integers(0, n_blocks)), block_len)
                else:
                    own = s61.adaptive_input(partner, x_out, Qs, mu, scale, Zsh,
                                             n_blocks, rng)
                Zp = s61.zscore(own, mu, scale)
                lam, V = s61.topk([Zh, Zp], k)
                Qs.append(V); Ys.append(V * lam[None, :] - (Zp.T @ Zp) @ V)
            q = t + 1
            if q in s61.QS:
                Q = np.concatenate(Qs, axis=1); Y = np.concatenate(Ys, axis=1)
                # the recovered operator Y Q^+ has the same column space as Y (Y = C_H Q
                # and Q has full row rank), so the membership score is the projection
                # norm onto span(Y) -- exactly s61.metrics' definition, and cheap
                # (svd of an M x q matrix, never M x M).
                auc_q, tpr_q, sgn_q = s61._membership_auc_from_operator(
                    Y, honest, nonmem, mu, scale, y, cal)
                per_q[str(q)] = {"auc": auc_q, "tpr05": tpr_q,
                                 "gap_vs_joint": auc_q - joint_auc,
                                 "gap_vs_joint_tpr05": tpr_q - joint_tpr,
                                 "gap_vs_tpr_oracle": tpr_q - tpr_oracle}
        out["strategies"][strategy] = per_q
    return out


def _boot_ci(diff):
    rng = np.random.default_rng(7)
    boot = diff[rng.integers(0, len(diff), size=(10000, len(diff)))].mean(1)
    return [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]


def run(chrom=22, seeds=5):
    out = {"chrom": chrom, "seeds": seeds, "qs": list(s61.QS),
           "ci_definition": "10k percentile bootstrap over the 5 repeated splits of ONE "
                            "fixed panel; conditional split uncertainty, not "
                            "independent-subject generalization",
           "cohorts": {}}
    for name, (G, groups) in cohorts_for(chrom).items():
        rows = [one_split(G, groups, s) for s in range(seeds)]
        blk = {"strategies": {}, "oracle_key_counts": {}, "per_split": rows}
        for key in ("joint_passive_auc", "joint_passive_tpr05", "tpr_oracle_bar",
                    "fused_z_auc", "fused_z_tpr05", "fused_cal_auc", "fused_cal_tpr05",
                    "pooled_center_scale_auc", "af_cond_auc", "spec_cond_auc"):
            blk[key + "_mean"] = float(np.mean([r[key] for r in rows]))
        for r in rows:
            blk["oracle_key_counts"][r["oracle_key"]] = \
                blk["oracle_key_counts"].get(r["oracle_key"], 0) + 1
            blk["tpr_oracle_key_counts"] = blk.get("tpr_oracle_key_counts", {})
            blk["tpr_oracle_key_counts"][r["tpr_oracle_key"]] = \
                blk["tpr_oracle_key_counts"].get(r["tpr_oracle_key"], 0) + 1
        for kk in rows[0]["joint_candidates_auc"]:
            blk.setdefault("joint_candidates_auc_mean", {})[kk] = float(
                np.mean([r["joint_candidates_auc"][kk] for r in rows]))
            blk.setdefault("joint_candidates_tpr05_mean", {})[kk] = float(
                np.mean([r["joint_candidates_tpr05"][kk] for r in rows]))
        for st in ("passive", "random", "adaptive"):
            blk["strategies"][st] = {}
            for q in s61.QS:
                auc = np.array([r["strategies"][st][str(q)]["auc"] for r in rows])
                tpr = np.array([r["strategies"][st][str(q)]["tpr05"] for r in rows])
                gap = np.array([r["strategies"][st][str(q)]["gap_vs_joint"] for r in rows])
                gap_t = np.array([r["strategies"][st][str(q)]["gap_vs_joint_tpr05"]
                                  for r in rows])
                gap_to = np.array([r["strategies"][st][str(q)]["gap_vs_tpr_oracle"]
                                   for r in rows])
                blk["strategies"][st][str(q)] = {
                    "auc": float(auc.mean()), "tpr05": float(tpr.mean()),
                    "gap_vs_joint": float(gap.mean()),
                    "gap_vs_joint_ci95": _boot_ci(gap),
                    "gap_vs_joint_tpr05": float(gap_t.mean()),
                    "gap_vs_joint_tpr05_ci95": _boot_ci(gap_t),
                    "pos_splits": int((gap > 0).sum()),
                    "pos_splits_tpr05": int((gap_t > 0).sum()),
                    "gap_vs_tpr_oracle": float(gap_to.mean()),
                    "gap_vs_tpr_oracle_ci95": _boot_ci(gap_to),
                    "pos_splits_tpr_oracle": int((gap_to > 0).sum())}
        out["cohorts"][name] = blk
        print("### %s  auc-oracle bar=%s tpr-oracle bar=%s  auc %.3f tpr %.3f  "
              "fused_z %.3f  fused_cal %.3f  pooled %.3f  af_cond %.3f  spec_cond %.3f" % (
                  name, blk["oracle_key_counts"], blk["tpr_oracle_key_counts"],
                  blk["joint_passive_auc_mean"], blk["joint_passive_tpr05_mean"],
                  blk["fused_z_auc_mean"], blk["fused_cal_auc_mean"],
                  blk["pooled_center_scale_auc_mean"], blk["af_cond_auc_mean"],
                  blk["spec_cond_auc_mean"]), flush=True)
        for st in ("passive", "random", "adaptive"):
            print("  %-8s " % st + " | ".join(
                "q%d auc%.3f gapA%+.4f [%+.4f,%+.4f] pos%d gapT%+.4f [%+.4f,%+.4f] posT%d" % (
                    q, blk["strategies"][st][str(q)]["auc"],
                    blk["strategies"][st][str(q)]["gap_vs_joint"],
                    blk["strategies"][st][str(q)]["gap_vs_joint_ci95"][0],
                    blk["strategies"][st][str(q)]["gap_vs_joint_ci95"][1],
                    blk["strategies"][st][str(q)]["pos_splits"],
                    blk["strategies"][st][str(q)]["gap_vs_tpr_oracle"],
                    blk["strategies"][st][str(q)]["gap_vs_tpr_oracle_ci95"][0],
                    blk["strategies"][st][str(q)]["gap_vs_tpr_oracle_ci95"][1],
                    blk["strategies"][st][str(q)]["pos_splits_tpr_oracle"]) for q in s61.QS),
                flush=True)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s67_joint_passive_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
