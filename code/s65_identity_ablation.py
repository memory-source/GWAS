"""S65 - identity ablation for the selection-input spectral tomography (deep round).

The round-2 review ran an ablation showing the tomography identity is load-bearing:
stacking the raw released eigenvectors V_t (without forming Y_t = V_t Lambda_t - K_t V_t)
does NOT beat the strong passive baseline, while the identity-based stack does. That
result was produced in the reviewer's scratch code; this script makes it a
machine-readable artifact.

For each cohort the script runs the SAME random block-swap sequence as s61.cell, keeps
the released V_t at each release, and scores membership two ways at q = 1/2/4/8:
  identity   -- the s61 stack, Y_t = V_t Lambda_t - K_t V_t (control path)
  raw_V      -- the union span of the released V_t alone
and reports the paired gap against the strongest passive baseline in each case.

Output: results/raw/s65_identity_ablation_chr<C>.json
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


def one_split(G, groups, seed, k=s61.K_TOP, block_len=s61.BLOCK_LEN, caps=None):
    """Re-run the random strategy's releases, keeping the released V_t, and score the
    identity stack and the raw-V stack at each q. Returns (gaps_identity, gaps_raw,
    strong)."""
    roles = s61.split_roles(groups, caps or s61.CAPS, seed)
    honest, partner = G[roles["honest"]], G[roles["partner"]]
    nonmem, shadow = G[roles["nonmember"]], G[roles["shadow"]]
    mu, scale = s61.joint_baseline(honest, partner)
    Zh = s61.zscore(honest, mu, scale)
    Zsh = s61.zscore(shadow, mu, scale)
    n_blocks = G.shape[1] // block_len
    y = np.concatenate([np.ones(len(honest)), np.zeros(len(nonmem))])
    cal = np.random.default_rng(seed + 31).random(len(y)) < 0.5
    grp_rows = np.concatenate([np.asarray(groups)[roles["honest"]],
                               np.asarray(groups)[roles["nonmember"]]])

    # the same passive first release defines the shared output 1 and the baseline
    Zp0 = s61.zscore(partner, mu, scale)
    lam0, V0 = s61.topk([Zh, Zp0], k)
    Y0 = V0 * lam0[None, :] - (Zp0.T @ Zp0) @ V0
    _, _, sgn_spec = s61._membership_auc_from_operator(Y0, honest, nonmem, mu, scale, y, cal)
    Zk = Zh[np.asarray(cal, bool)[:len(honest)]]
    Y0c = Y0 - Zk.T @ (Zk @ V0)
    spec_cond, _, _ = s61._membership_auc_from_operator(
        Y0c, honest, nonmem, mu, scale, y, cal, sign=sgn_spec)
    af_cond = s61.af_baselines(honest, nonmem, G[roles["reference"]], partner, shadow,
                               np.asarray(groups)[roles["honest"]], y, cal, seed)[2]
    strong = max(spec_cond, af_cond)

    rng = np.random.default_rng(seed + 101)
    Vs, Ys = [V0], [Y0]
    gap_id, gap_raw = {}, {}
    # q = 1 is the shared passive first release (both stacks see the same V_1)
    for tag, mats in (("id", Ys), ("raw", Vs)):
        auc, _, _ = s61._membership_auc_from_operator(mats[0], honest, nonmem, mu, scale, y, cal)
        (gap_id if tag == "id" else gap_raw)["1"] = auc - strong
    for t in range(1, max(s61.QS)):
        own, _ = s61.swap_block(partner, rng.choice(partner.shape[0], 2, replace=False),
                                int(rng.integers(0, n_blocks)), block_len)
        Zp = s61.zscore(own, mu, scale)
        lam, V = s61.topk([Zh, Zp], k)
        Vs.append(V); Ys.append(V * lam[None, :] - (Zp.T @ Zp) @ V)
        q = t + 1
        if q in s61.QS:
            Qid = np.concatenate(Ys[:q], axis=1)
            auc_id, _, _ = s61._membership_auc_from_operator(Qid, honest, nonmem, mu, scale, y, cal)
            Qraw = np.concatenate(Vs[:q], axis=1)
            auc_raw, _, _ = s61._membership_auc_from_operator(Qraw, honest, nonmem, mu, scale, y, cal)
            gap_id[str(q)] = auc_id - strong
            gap_raw[str(q)] = auc_raw - strong
    return gap_id, gap_raw, strong, grp_rows


def run(chrom=22, seeds=5):
    out = {"chrom": chrom, "seeds": seeds, "qs": list(s61.QS), "cohorts": {}}
    for name, (G, groups) in cohorts_for(chrom).items():
        rows = [one_split(G, groups, s) for s in range(seeds)]
        blk = {"strategies": {}}
        for tag, idx in (("identity", 0), ("raw_V", 1)):
            blk["strategies"][tag] = {}
            for q in s61.QS:
                g = np.array([r[idx][str(q)] for r in rows])
                blk["strategies"][tag][str(q)] = {
                    "gap_vs_strongest_passive": float(g.mean()),
                    "gap_vs_strongest_passive_ci95": s61._boot_ci(g),
                    "pos_splits": int((g > 0).sum())}
        blk["strongest_passive_mean"] = float(np.mean([r[2] for r in rows]))
        out["cohorts"][name] = blk
        print("### %s  strongest passive %.3f" % (name, blk["strongest_passive_mean"]),
              flush=True)
        for tag in ("identity", "raw_V"):
            print("  %-9s " % tag + " | ".join(
                "q%d gap%+.4f" % (q, blk["strategies"][tag][str(q)]["gap_vs_strongest_passive"])
                for q in s61.QS), flush=True)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s65_identity_ablation_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
