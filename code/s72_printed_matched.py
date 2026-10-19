"""S72 - matched printed-output comparison (round 22, closing the J9 gap).

The round-22 review (J9) asked for the decisive matched comparison: measure the
accumulated active attack against the passive adversary built from the SAME PRINTED
PLINK outputs, not float64 memory. This script runs, per cohort x split seed, a
sequence of q = 8 joint re-fits through real PLINK 2.0 (`--pca <k> allele-wts approx`),
reads the printed allele weights and eigenvalues back, forms the attack operand
`Y_t = V_t Λ_t − K_t V_t` from the READ-BACK V/Λ (with the partner's own scatter K_t
known exactly), and compares the q8 stack against the four passive constructions
(af_cond, spec_cond, fused_z, fused_cal) built from the same read-back outputs.

af_cond uses the released pooled AF (the plan grants it); spec_cond uses the read-back
single release's first output conditioned by the calibration members' scatter.

Output: results/raw/s72_printed_matched_chr<C>.json
"""
import json
import os
import subprocess
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
import s51_real_lowstruct as s51
import s61_spectral_tomography as s61
import s63_plink_export as s63
from common import allele_freq, af_lrt_attack, auroc, tpr_at_fpr
from real_multi import SUPERPOP

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "raw")
WORK = os.path.join(ROOT, "scratch", "printed")
K = 10
QS = (1, 2, 4, 8)


def _cohorts(chrom, ld_window=1500, max_snps=20000):
    Ge, _, pos = s48.load_for(("EUR",), chrom, ld_window, max_snps)
    G3, _, sup3, pos3 = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    return {"within:EUR": (Ge, np.array(["EUR"] * len(Ge)), pos),
            "super:3pop": (G3, sup3, pos3)}


def _z(x):
    sd = x.std()
    return (x - x.mean()) / (sd if sd > 0 else 1.0)


def _pnorm(R, honest, nonmem, mu, scale):
    U, s, _ = np.linalg.svd(R, full_matrices=False)
    keep = s > 1e-9 * (s[0] if s.size else 1.0)
    X = s61.zscore(np.vstack([honest, nonmem]), mu, scale)
    return np.linalg.norm(X @ U[:, keep], axis=1)


def _spec_cond(Y, Zk, V, honest, nonmem, mu, scale, M):
    """Conditioned single release: subtract the known calibration members' projection
    onto the released directions V (Zk @ V, in PLINK's scatter/M normalization), then
    score by the projection norm onto the span of the remainder."""
    Yc = Y - Zk.T @ (Zk @ V) / M
    return _pnorm(Yc, honest, nonmem, mu, scale)


def one_split(G, groups, seed, pos, chrom, cohort, block_len=s61.BLOCK_LEN, caps=None):
    n_rows = G.shape[0]
    joint_ids = [f"I{i}" for i in range(n_rows)]              # names do not affect PCA
    roles = s61.split_roles(groups, caps or s61.CAPS, seed)
    honest, partner = G[roles["honest"]], G[roles["partner"]]
    nonmem, ref, shadow = G[roles["nonmember"]], G[roles["reference"]], G[roles["shadow"]]
    mu, scale = s61.joint_baseline(honest, partner)
    Zh = s61.zscore(honest, mu, scale)
    n_blocks = G.shape[1] // block_len
    y = np.concatenate([np.ones(len(honest)), np.zeros(len(nonmem))])
    cal = np.random.default_rng(seed + 31).random(len(y)) < 0.5
    h = np.asarray(cal, bool)
    honest_groups = np.asarray(groups)[roles["honest"]]

    rng = np.random.default_rng(seed + 101)
    Ys, Vs, Zps64, Vs64, lams64 = [], [], [], [], []
    V0, lam0 = None, None
    M = G.shape[1]
    for t in range(max(QS)):
        if t == 0:
            own = partner
        else:
            own, _ = s61.swap_block(partner, rng.choice(partner.shape[0], 2, replace=False),
                                    int(rng.integers(0, n_blocks)), block_len)
        jp = np.vstack([honest, own])
        vcf = os.path.join(WORK, f"s72_c{chrom}_{cohort}_{seed % 10000}_{t}.vcf")
        s63.write_vcf(vcf, jp, pos, [f"H{i}" for i in range(len(honest))] +
                      [f"P{i}" for i in range(len(own))])
        pref = os.path.join(WORK, f"s72_c{chrom}_{cohort}_{seed % 10000}_{t}")
        cmd = [s63.PLINK, "--vcf", vcf, "--allow-no-sex", "--pca", str(K), "allele-wts",
               "approx", "--seed", str(1000 + t), "--threads", "1", "--out", pref]
        pr = subprocess.run(cmd, capture_output=True, text=True)
        if pr.returncode != 0:
            raise RuntimeError(f"plink2 rc={pr.returncode}\n{pr.stdout[-500:]}")
        lam, V, _ids = s63.read_pca(pref)
        Zp = s61.zscore(own, mu, scale)
        K_mat = Zp.T @ Zp / M                           # PLINK eigenvalues are scatter/M
        Ys.append(V * lam[None, :] - K_mat @ V)
        Vs.append(V)
        Zps64.append(Zp)
        lam64, V64 = s61.topk([Zh, Zp], K)              # float64 reference in the same scaling
        Vs64.append(V64)
        lams64.append(lam64)
        if t == 0:
            V0, lam0 = V, lam
    Yq = np.concatenate(Ys, axis=1)

    # ---- passive channels from the same read-back outputs ----
    # The passive comparator is built under the SAME information model as s67's spec_cond:
    # it sees the same released V/Λ and the attack operand Y_t = V_t Λ_t − K_t V_t (the
    # conditioned single release subtracts the known calibration members' projection
    # Zk@V_t from Y_t), and fuses the per-release scalar channels over q. This is a
    # late-fusion diagnostic over the specified constructions, not a general statement
    # about arbitrary processing of the release transcript.
    af_s = allele_freq(honest)
    bg = allele_freq(np.vstack([ref, partner, shadow]))
    raw_pooled = np.concatenate([af_lrt_attack(honest, af_s, bg),
                                 af_lrt_attack(nonmem, af_s, bg)])
    quota = min(int((h[:len(honest)] & (honest_groups == g)).sum())
                for g in np.unique(honest_groups))
    bal = np.concatenate([np.where(h[:len(honest)] & (honest_groups == g))[0][:quota]
                          for g in np.unique(honest_groups)])
    num = 2.0 * len(honest) * af_s - honest[bal].sum(0).astype(float)
    af_c = np.clip(num / (2.0 * (len(honest) - len(bal))), 1e-6, 1 - 1e-6)
    raw_cond = np.concatenate([af_lrt_attack(honest, af_c, bg),
                               af_lrt_attack(nonmem, af_c, bg)])
    sgn_af = 1.0 if auroc(raw_pooled[h], y[h]) >= 0.5 else -1.0
    s_af_cond = sgn_af * raw_cond
    s_pooled = sgn_af * raw_pooled

    def heldout(s):
        return float(auroc(s[~h], y[~h])), float(tpr_at_fpr(s[~h], y[~h], 0.05))

    from sklearn.linear_model import LogisticRegression

    def oracle_bar(spec_cond_channels, spec_uncond_channels):
        """Best held-out AUROC over the four passive constructions of one output budget:
        af_cond, the spec channel(s) fused if several, fused_z (label-free), and
        fused_cal (logistic on the calibration half). Mirroring s67, the fused_cal
        weights are fitted on the UNCONDITIONED calibration features (s_pooled, the
        unconditioned spec) and applied to the CONDITIONED channels, so the conditioning
        that removes the calibration members does not also destroy the fitted weights."""
        spec_c = np.sum([_z(c) for c in spec_cond_channels], axis=0)
        spec_u = np.sum([_z(c) for c in spec_uncond_channels], axis=0)
        cands = {"af_cond": s_af_cond, "spec_cond": spec_c,
                 "fused_z": _z(s_af_cond) + _z(spec_c)}
        Z_uncond = np.column_stack([_z(s_pooled), _z(spec_u)])
        w = LogisticRegression(max_iter=1000).fit(Z_uncond[h], y[h]).coef_[0]
        cands["fused_cal"] = w[0] * _z(s_af_cond) + w[1] * _z(spec_c)
        ho = {k: heldout(v) for k, v in cands.items()}
        key = max(ho, key=lambda k: ho[k][0])
        return key, ho[key][0], ho[key][1], {k: ho[k][0] for k in ho}

    # single-release passive bar (t = 0 only) and a late-fusion bar over the q outputs
    Zk = Zh[h[:len(honest)]]
    spec_single = _spec_cond(Ys[0], Zk, V0, honest, nonmem, mu, scale, G.shape[1])
    spec_uncond_single = _pnorm(Ys[0], honest, nonmem, mu, scale)
    spec_all = [_spec_cond(Yt, Zk, Vt, honest, nonmem, mu, scale, G.shape[1])
                for Yt, Vt in zip(Ys, Vs)]
    spec_uncond_all = [_pnorm(Yt, honest, nonmem, mu, scale) for Yt in Ys]
    key1, auc1, tpr1, cands1 = oracle_bar([spec_single], [spec_uncond_single])
    keyq, aucq, tprq, candsq = oracle_bar(spec_all, spec_uncond_all)

    # q8 stack from the read-back outputs, oriented on the calibration half
    auc_q8, tpr_q8, _ = s61._membership_auc_from_operator(Yq, honest, nonmem, mu, scale,
                                                          y, cal)
    return {"seed": seed,
            "oracle_key_single": key1, "oracle_auc_single": auc1,
            "oracle_tpr05_single": tpr1,
            "oracle_key_q": keyq, "oracle_auc_q": aucq, "oracle_tpr05_q": tprq,
            "q8_auc": auc_q8, "q8_tpr05": tpr_q8,
            "gap_q8_vs_oracle_single": auc_q8 - auc1,
            "gap_tpr05_vs_oracle_single": tpr_q8 - tpr1,
            "gap_q8_vs_oracle_q": auc_q8 - aucq,
            "gap_tpr05_vs_oracle_q": tpr_q8 - tprq,
            "candidates_auc_single": cands1, "candidates_auc_q": candsq}


def _boot_ci(diff):
    rng = np.random.default_rng(7)
    boot = diff[rng.integers(0, len(diff), size=(10000, len(diff)))].mean(1)
    return [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]


def run(chrom=22, seeds=5):
    out = {"chrom": chrom, "seeds": seeds, "k": K, "qs": list(QS),
           "solver": "PLINK 2.0 --pca allele-wts approx", "cohorts": {}}
    for name, (G, groups, pos) in _cohorts(chrom).items():
        tag = name.split(":")[-1]
        rows = [one_split(G, groups, s, pos, chrom, tag) for s in range(seeds)]
        gap1 = np.array([r["gap_q8_vs_oracle_single"] for r in rows])
        gap1_t = np.array([r["gap_tpr05_vs_oracle_single"] for r in rows])
        gapq = np.array([r["gap_q8_vs_oracle_q"] for r in rows])
        gapq_t = np.array([r["gap_tpr05_vs_oracle_q"] for r in rows])
        out["cohorts"][name] = {
            "per_split": rows,
            "oracle_auc_single_mean": float(np.mean([r["oracle_auc_single"] for r in rows])),
            "oracle_auc_q_mean": float(np.mean([r["oracle_auc_q"] for r in rows])),
            "q8_auc_mean": float(np.mean([r["q8_auc"] for r in rows])),
            "gap_q8_vs_oracle_single": float(gap1.mean()),
            "gap_q8_vs_oracle_single_ci95": _boot_ci(gap1),
            "gap_tpr05_vs_oracle_single": float(gap1_t.mean()),
            "gap_tpr05_vs_oracle_single_ci95": _boot_ci(gap1_t),
            "pos_splits_single": int((gap1 > 0).sum()),
            "gap_q8_vs_oracle_q": float(gapq.mean()),
            "gap_q8_vs_oracle_q_ci95": _boot_ci(gapq),
            "gap_tpr05_vs_oracle_q": float(gapq_t.mean()),
            "gap_tpr05_vs_oracle_q_ci95": _boot_ci(gapq_t),
            "pos_splits_q": int((gapq > 0).sum())}
        b = out["cohorts"][name]
        print("### %s  oracle1 %.3f oracle_q %.3f (q8 %.3f)  gap1 %+.4f pos %d  "
              "gapq %+.4f pos %d" % (
                  name, b["oracle_auc_single_mean"], b["oracle_auc_q_mean"],
                  b["q8_auc_mean"], b["gap_q8_vs_oracle_single"], b["pos_splits_single"],
                  b["gap_q8_vs_oracle_q"], b["pos_splits_q"]), flush=True)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s72_printed_matched_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
