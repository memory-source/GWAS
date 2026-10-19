"""S70 - known-member-fraction sweep (round 22).

The calibrated passive adversary and the conditioned attack both use the KNOWN
calibration members. The GPT analysis asks how the increment behaves as that knowledge
shrinks (0% / 5% / 10% / 50%). The known knowledge at fraction f is: f * N honest
members (drawn from the calibration half) AND an equal number of NON-members (drawn
from the calibration half's nonmembers), so the label set is two-class and can fix an
orientation. At 0% no member label exists and every orientation uses the mechanism
sign (+1: a member's scatter entered C_H, so it projects further along the released
range) -- an a-priori direction from the identity, not label leakage.

For each cohort x split seed x fraction:
  - re-form af_cond (composition-balanced removal of the known members from the
    released AF) and spec_cond (the single release with the known members' scatter
    removed);
  - orient on the known label set (or the mechanism sign at 0%), evaluate held-out;
  - rebuild the fused_z passive adversary and the q=8 random stack; report the paired
    gap with bootstrap CI and the positive split count.

Output: results/raw/s70_known_fraction_chr<C>.json
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
import s67_joint_passive as s67
from common import allele_freq, af_lrt_attack, auroc, tpr_at_fpr
from real_multi import SUPERPOP

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
FRACTIONS = (0.0, 0.05, 0.10, 0.50)


def cohorts_for(chrom, ld_window=1500, max_snps=20000):
    Ge, _, _ = s48.load_for(("EUR",), chrom, ld_window, max_snps)
    G3, _, sup3, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    G5, _, sup5, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    return {"within:EUR": (Ge, np.array(["EUR"] * len(Ge))),
            "super:3pop": (G3, sup3),
            "within:EAS": (G5[sup5 == "EAS"], sup5[sup5 == "EAS"])}


def _pnorm(R, honest, nonmem, mu, scale):
    U, s, _ = np.linalg.svd(R, full_matrices=False)
    keep = s > 1e-9 * (s[0] if s.size else 1.0)
    X = s61.zscore(np.vstack([honest, nonmem]), mu, scale)
    return np.linalg.norm(X @ U[:, keep], axis=1)


def one_split(G, groups, seed, k=s61.K_TOP, block_len=s61.BLOCK_LEN, caps=None, qs=(8,)):
    roles = s61.split_roles(groups, caps or s61.CAPS, seed)
    honest, partner = G[roles["honest"]], G[roles["partner"]]
    nonmem, ref, shadow = G[roles["nonmember"]], G[roles["reference"]], G[roles["shadow"]]
    mu, scale = s61.joint_baseline(honest, partner)
    Zh = s61.zscore(honest, mu, scale)
    n_blocks = G.shape[1] // block_len
    N = len(honest)
    y = np.concatenate([np.ones(N), np.zeros(len(nonmem))])
    cal = np.random.default_rng(seed + 31).random(len(y)) < 0.5
    h = np.asarray(cal, bool)
    honest_groups = np.asarray(groups)[roles["honest"]]

    # random release sequence, reusing the s61 swap generator, keeping V_t and Y_t
    Zp0 = s61.zscore(partner, mu, scale)
    lam0, V0 = s61.topk([Zh, Zp0], k)
    Y0 = V0 * lam0[None, :] - (Zp0.T @ Zp0) @ V0
    rng = np.random.default_rng(seed + 101)
    Ys = [Y0]
    for t in range(1, max(qs)):
        own, _ = s61.swap_block(partner, rng.choice(partner.shape[0], 2, replace=False),
                                int(rng.integers(0, n_blocks)), block_len)
        Zp = s61.zscore(own, mu, scale)
        lam, V = s61.topk([Zh, Zp], k)
        Ys.append(V * lam[None, :] - (Zp.T @ Zp) @ V)
    Yq = np.concatenate(Ys[:max(qs)], axis=1)

    af_s = allele_freq(honest)
    bg_pooled = allele_freq(np.vstack([ref, partner, shadow]))
    raw_pooled = np.concatenate([af_lrt_attack(honest, af_s, bg_pooled),
                                 af_lrt_attack(nonmem, af_s, bg_pooled)])
    d_spec_single = _pnorm(Y0, honest, nonmem, mu, scale)
    d_q = _pnorm(Yq, honest, nonmem, mu, scale)

    # calibration-half member pool (label-access budget for the CONDITIONING only)
    cal_mem = np.where(h & (y == 1))[0]

    out = {"seed": seed, "fractions": {}}
    for frac in FRACTIONS:
        # the fraction is a share of the HONEST COHORT (N members); the drawn members
        # are restricted to the calibration half (the labelled part)
        n_mem = int(round(frac * N))
        rngk = np.random.default_rng(seed + 200)
        mem_known = rngk.choice(cal_mem, size=min(n_mem, len(cal_mem)), replace=False)
        known_members_h = np.zeros(N, bool)
        if len(mem_known) > 0:
            known_members_h[mem_known] = True
        label_free = (len(mem_known) == 0)

        # Orientation is the MECHANISM sign (+1) for both channels at EVERY fraction:
        # a study member's genotype matches the released study AF better (AF-LRT score
        # higher) and projects further along the released range (its scatter entered
        # C_H). This is an a-priori direction from the identity, so it needs NO member
        # labels. The known-member fraction therefore controls ONLY the conditioning
        # (which members are subtracted), never the orientation.
        def ev(v):
            sc = v[~h]
            return float(auroc(sc, y[~h])), float(tpr_at_fpr(sc, y[~h], 0.05))

        # conditioned AF: remove an equal number of known members from EACH group
        # (composition-balanced: the quota is the smallest known count across groups)
        if known_members_h.sum() > 0:
            gh = honest_groups
            counts = [int((known_members_h & (gh == g)).sum()) for g in np.unique(gh)]
            quota = min(counts) if counts else 0
            if quota > 0:
                bal = np.concatenate([np.where(known_members_h & (gh == g))[0][:quota]
                                      for g in np.unique(gh)])
            else:
                bal = np.array([], int)
            if len(bal) == 0:
                bal = np.where(known_members_h)[0]
            num = 2.0 * N * af_s - honest[bal].sum(0).astype(float)
            af_c = np.clip(num / (2.0 * (N - len(bal))), 1e-6, 1 - 1e-6)
        else:
            bal = np.array([], int)
            af_c = af_s
        raw_c = np.concatenate([af_lrt_attack(honest, af_c, bg_pooled),
                                af_lrt_attack(nonmem, af_c, bg_pooled)])

        # conditioned single release
        if known_members_h.sum() > 0:
            Zk = Zh[known_members_h]
            Y0c = Y0 - Zk.T @ (Zk @ V0)
        else:
            Y0c = Y0
        d_spec_c = _pnorm(Y0c, honest, nonmem, mu, scale)

        # the passive bar is the LABEL-FREE z-sum fusion (reachable at every fraction
        # including 0% known, since the orientation is the mechanism sign). No
        # label-fitted fusion is reported here: at a given known-member budget the
        # attacker holds only that many labels, so a fusion fitted on all calibration
        # labels would exceed the budget.
        fused = s67._z(raw_c) + s67._z(d_spec_c)
        joint_auc, joint_tpr = ev(fused)
        q8_auc, q8_tpr = ev(d_q)
        out["fractions"][str(frac)] = {
            "n_known_members": int(len(mem_known)), "n_balanced_removed": int(len(bal)),
            "label_free": bool(label_free),
            "joint_passive_auc": joint_auc, "joint_passive_tpr05": joint_tpr,
            "q8_auc": q8_auc, "q8_tpr05": q8_tpr,
            "q8_auc_gap_vs_tag": "label_free_z_fusion",
            "gap_q8_vs_joint": q8_auc - joint_auc,
            "gap_q8_tpr05": q8_tpr - joint_tpr}
    return out


def _boot_ci(diff):
    rng = np.random.default_rng(7)
    boot = diff[rng.integers(0, len(diff), size=(10000, len(diff)))].mean(1)
    return [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]


def run(chrom=22, seeds=5):
    out = {"chrom": chrom, "seeds": seeds, "fractions": list(FRACTIONS),
           "ci_definition": "10k percentile bootstrap over the 5 repeated splits of ONE "
                            "fixed panel; conditional split uncertainty, not "
                            "independent-subject generalization",
           "cohorts": {}}
    for name, (G, groups) in cohorts_for(chrom).items():
        rows = [one_split(G, groups, s) for s in range(seeds)]
        blk = {"per_split": rows,
               "n_known_members_per_split": [r["fractions"][str(FRACTIONS[3])]["n_known_members"]
                                             for r in rows]}
        for frac in FRACTIONS:
            fr = str(frac)
            auc = np.array([r["fractions"][fr]["q8_auc"] for r in rows])
            gap = np.array([r["fractions"][fr]["gap_q8_vs_joint"] for r in rows])
            jp = np.array([r["fractions"][fr]["joint_passive_auc"] for r in rows])
            fcs = np.array([r["fractions"][fr].get("fused_cal_auc", np.nan) for r in rows], float)
            blk[fr] = {"joint_passive_auc": float(jp.mean()),
                       "q8_auc": float(auc.mean()),
                       "gap_q8_vs_joint": float(gap.mean()),
                       "gap_q8_vs_joint_ci95": _boot_ci(gap),
                       "pos_splits": int((gap > 0).sum()),
                       "n_known_members_mean": float(np.mean(
                           [r["fractions"][fr]["n_known_members"] for r in rows])),
                       "n_balanced_removed_mean": float(np.mean(
                           [r["fractions"][fr]["n_balanced_removed"] for r in rows]))}
        out["cohorts"][name] = blk
        print("### %s  n_known@50%% %s" % (name, blk["n_known_members_per_split"]), flush=True)
        for frac in FRACTIONS:
            b = blk[str(frac)]
            print("   known %.2f  joint %.3f  q8 %.3f  gap %+.4f [%+.4f,%+.4f] pos %d" % (
                frac, b["joint_passive_auc"], b["q8_auc"], b["gap_q8_vs_joint"],
                b["gap_q8_vs_joint_ci95"][0], b["gap_q8_vs_joint_ci95"][1],
                b["pos_splits"]), flush=True)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s70_known_fraction_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
