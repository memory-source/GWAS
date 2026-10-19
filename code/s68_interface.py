"""S68 - real interface manifest + approximate-solver test (round 22).

Two questions that gate the whole route:

1. INTERFACE MANIFEST. The attack assumes the attacker can obtain, per re-fit, the
   joint top-k loadings, the eigenvalues, the public center/scale, and the right to
   re-run on the same honest cohort after modifying its own input, at the printed
   precision. This script states each item as {obtainable, demonstrated_by, note}
   and demonstrates the loadings/eigenvalues/precision items with a real PLINK 2.0
   `--pca <k> allele-wts` read-back.

2. APPROXIMATE SOLVER. PLINK's `--pca approx` replaces the exact deterministic
   eigen-decomposition with the randomized Galinsky-2016 algorithm. The
   exact-equation attack may not survive the solver residual. This script runs the
   joint PCA with `--pca <k> allele-wts approx` on the same joint VCF, reads the
   result back, and computes (i) the identity residual, (ii) the q8 stacked attack
   vs the conditioned single release, and (iii) the top-k span movement between two
   DIFFERENT solver seeds on the SAME input (the solver-randomness control).

Output: results/raw/s68_interface_chr<C>.json
"""
import json
import os
import subprocess
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
import s61_spectral_tomography as s61
import s63_plink_export as s63
from real_multi import load_panel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "raw")
PLINK = os.path.join(ROOT, "scratch", "plink2bin", "plink2")
WORK = os.path.join(ROOT, "scratch", "interface")
K = 10
BLOCK_LEN = s61.BLOCK_LEN


def run_plink_pca(vcf, out_prefix, k=K, approx=False, seed=None):
    cmd = [PLINK, "--vcf", vcf, "--allow-no-sex", "--pca", str(k), "allele-wts"]
    if approx:
        cmd.append("approx")
    if seed is not None:
        cmd += ["--seed", str(seed), "--threads", "1"]
    cmd += ["--out", out_prefix]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"plink2 failed rc={p.returncode}\n{p.stdout[-800:]}\n{p.stderr[-800:]}")


def read_pca(out_prefix, k=K):
    with open(out_prefix + ".eigenval") as f:
        lam = np.array([float(x) for x in f.read().split()])
    ids, rows = [], []
    with open(out_prefix + ".eigenvec.allele") as f:
        f.readline()
        for line in f:
            p = line.split()
            ids.append(p[1])
            rows.append([float(x) for x in p[5:]])
    rows = np.array(rows)
    return lam[:k], rows[1::2], ids[1::2]


def subspace_distance(A, B):
    """Principal angles between two column spaces; returns max sin(principal angle)."""
    Qa = np.linalg.qr(A)[0]
    Qb = np.linalg.qr(B)[0]
    s = np.linalg.svd(Qa.T @ Qb, compute_uv=False)
    s = np.clip(s, -1.0, 1.0)
    return float(np.sqrt(max(0.0, 1.0 - s.min() ** 2)))


def run(chrom=22, seeds=3):
    G, lab, pos = s48.load_for(("EUR",), chrom, 1500, 20000)
    groups = np.array(["EUR"] * len(G))
    samples, sups = load_panel()
    sample_ids = [s for s, k in zip(samples, sups) if k == "EUR"]
    os.makedirs(WORK, exist_ok=True)
    n_blocks = G.shape[1] // BLOCK_LEN
    M = G.shape[1]

    manifest = {
        "joint_topk_loadings": {
            "obtainable": True,
            "demonstrated_by": "PLINK 2.0 --pca <k> allele-wts -> .eigenvec.allele (ALT rows)",
            "note": "the joint PCA is on (honest, partner); the partner receives the "
                    "returned allele weights"},
        "eigenvalues": {
            "obtainable": True,
            "demonstrated_by": "PLINK 2.0 --pca <k> allele-wts -> .eigenval",
            "note": "equals scatter-eigenvalues / marker count"},
        "public_center_scale": {
            "obtainable": True,
            "demonstrated_by": "released pooled AF (honest+partner mean)",
            "note": "deterministic function of the released AF; grants nothing beyond it"},
        "rerun_right": {
            "obtainable": True,
            "demonstrated_by": "joint VCF rewritten with the partner's modified input per release",
            "note": "the same honest rows stay in the joint VCF; only the partner's own "
                    "rows change"},
        "printed_precision": {
            "obtainable": True,
            "demonstrated_by": "the .eigenvec.allele decimals (6 significant digits)",
            "note": "measured as the identity residual below"},
        "honest_private_matrix": {
            "obtainable": False,
            "demonstrated_by": "n/a",
            "note": "the attack never reads C_H or the honest genotypes beyond the "
                    "released output and the public center/scale"},
    }

    per_seed = []
    for s in range(seeds):
        roles = s61.split_roles(groups, s61.CAPS, s)
        honest = G[roles["honest"]]; partner = G[roles["partner"]]
        nonmem = G[roles["nonmember"]]
        joint_ids = [sample_ids[i] for i in roles["honest"]] + \
                    [sample_ids[i] for i in roles["partner"]]
        p = np.vstack([honest, partner]).mean(0) / 2.0
        scale = np.sqrt(np.clip(2.0 * p * (1 - p), 1e-9, None))
        Zh = s63.std_panel(honest, p, scale)
        C_H = Zh.T @ Zh / M
        mu = 2.0 * p
        y = np.concatenate([np.ones(len(honest)), np.zeros(len(nonmem))])
        cal = np.random.default_rng(s + 31).random(len(y)) < 0.5

        # ---- exact solver: the s63 reference path ----
        rng = np.random.default_rng(s + 101)
        Qs, Ys = [], []
        resid_exact = 0.0
        order_ok = True
        for t in range(8):
            own = partner if t == 0 else s61.swap_block(
                partner, rng.choice(partner.shape[0], 2, replace=False),
                int(rng.integers(0, n_blocks)), BLOCK_LEN)[0]
            jp = np.vstack([honest, own])
            vcf = os.path.join(WORK, f"exact_s{s}_t{t}.vcf")
            s63.write_vcf(vcf, jp, pos, joint_ids)
            pref = os.path.join(WORK, f"exact_s{s}_t{t}")
            run_plink_pca(vcf, pref, approx=False)
            lam, V, ids_out = read_pca(pref)
            order_ok = order_ok and all(ids_out[j] == str(int(pos[j]))
                                        for j in range(len(ids_out)))
            Zp = s63.std_panel(own, p, scale)
            Kp = Zp.T @ Zp / M
            Y = V * lam[None, :] - Kp @ V
            resid_exact = max(resid_exact, float(np.max(np.abs(C_H @ V - Y)) / max(1e-12, np.max(np.abs(C_H @ V)))))
            Qs.append(V); Ys.append(Y)
        Qe = np.concatenate(Qs[:8], axis=1); Ye = np.concatenate(Ys[:8], axis=1)
        auc_exact, tpr_exact, _ = s61._membership_auc_from_operator(
            Ye @ np.linalg.pinv(Qe), honest, nonmem, mu, scale, y, cal)

        # ---- approximate solver on the same joint VCFs ----
        rng = np.random.default_rng(s + 101)
        Qa, Ya = [], []
        resid_approx = 0.0
        for t in range(8):
            own = partner if t == 0 else s61.swap_block(
                partner, rng.choice(partner.shape[0], 2, replace=False),
                int(rng.integers(0, n_blocks)), BLOCK_LEN)[0]
            jp = np.vstack([honest, own])
            vcf = os.path.join(WORK, f"approx_s{s}_t{t}.vcf")
            s63.write_vcf(vcf, jp, pos, joint_ids)
            pref = os.path.join(WORK, f"approx_s{s}_t{t}")
            run_plink_pca(vcf, pref, approx=True, seed=1000 + t)
            lam, V, ids_out = read_pca(pref)
            Zp = s63.std_panel(own, p, scale)
            Kp = Zp.T @ Zp / M
            Y = V * lam[None, :] - Kp @ V
            resid_approx = max(resid_approx, float(np.max(np.abs(C_H @ V - Y)) / max(1e-12, np.max(np.abs(C_H @ V)))))
            Qa.append(V); Ya.append(Y)
        Qap = np.concatenate(Qa[:8], axis=1); Yap = np.concatenate(Ya[:8], axis=1)
        auc_approx, tpr_approx, _ = s61._membership_auc_from_operator(
            Yap @ np.linalg.pinv(Qap), honest, nonmem, mu, scale, y, cal)

        # ---- solver-randomness control: two seeds on the SAME input ----
        own0 = partner
        jp = np.vstack([honest, own0])
        vcf = os.path.join(WORK, f"seedA_s{s}.vcf")
        s63.write_vcf(vcf, jp, pos, joint_ids)
        prefA = os.path.join(WORK, f"seedA_s{s}")
        prefB = os.path.join(WORK, f"seedB_s{s}")
        run_plink_pca(vcf, prefA, approx=True, seed=17)
        run_plink_pca(vcf, prefB, approx=True, seed=99991)
        _, VA, _ = read_pca(prefA)
        _, VB, _ = read_pca(prefB)
        span_move = subspace_distance(VA, VB)                    # sin of max angle
        # stacking the two seed spans on the same input: does it grow the member score?
        lamA, _, _ = read_pca(prefA)
        lamB, _, _ = read_pca(prefB)
        Zp0 = s63.std_panel(partner, p, scale); K0 = Zp0.T @ Zp0 / M
        YA = VA * lamA[None, :] - K0 @ VA
        YB = VB * lamB[None, :] - K0 @ VB
        _, tpr_stack, _ = s61._membership_auc_from_operator(
            np.concatenate([YA, YB], axis=1), honest, nonmem, mu, scale, y, cal)
        _, tpr_one, _ = s61._membership_auc_from_operator(YA, honest, nonmem, mu, scale, y, cal)

        per_seed.append({
            "seed": s, "marker_order_ok_approx": bool(
                all(ids_out[j] == str(int(pos[j])) for j in range(len(ids_out)))),
            "rel_identity_residual_exact": resid_exact,
            "rel_identity_residual_approx": resid_approx,
            "auc_q8_exact": auc_exact, "tpr05_q8_exact": tpr_exact,
            "auc_q8_approx": auc_approx, "tpr05_q8_approx": tpr_approx,
            "same_input_span_move_sin": span_move,
            "same_input_stack_tpr05": tpr_stack, "same_input_single_tpr05": tpr_one})
        print("seed %d  resid exact %.2e approx %.2e  auc exact %.3f approx %.3f  "
              "span_move %.3f  tpr stack %.3f vs single %.3f" % (
                  s, resid_exact, resid_approx, auc_exact, auc_approx, span_move,
                  tpr_stack, tpr_one), flush=True)

    res = {"chrom": chrom, "seeds": seeds, "k": K, "interface_manifest": manifest,
           "per_seed": per_seed,
           "resid_exact_max": float(np.max([r["rel_identity_residual_exact"] for r in per_seed])),
           "resid_approx_max": float(np.max([r["rel_identity_residual_approx"] for r in per_seed])),
           "auc_q8_exact_mean": float(np.mean([r["auc_q8_exact"] for r in per_seed])),
           "auc_q8_approx_mean": float(np.mean([r["auc_q8_approx"] for r in per_seed])),
           "span_move_mean": float(np.mean([r["same_input_span_move_sin"] for r in per_seed])),
           "same_input_stack_gain_mean": float(np.mean(
               [r["same_input_stack_tpr05"] - r["same_input_single_tpr05"] for r in per_seed]))}
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    ap.add_argument("--seeds", type=int, default=3)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom, a.seeds)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s68_interface_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
