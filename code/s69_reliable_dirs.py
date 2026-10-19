"""S69 - precision-reliable directions (round 22).

The round-19 stack used a FIXED relative rank cutoff (1e-9) on the singular values of
the released operand stack Y to decide which directions the attack uses. A direction
that is numerically new is not necessarily usable once the export is rounded. This
script derives a reliability floor from the MEASURED precision-induced perturbation of
Y and re-measures coverage and membership on the reliable set only.

The membership attack (s61.metrics) scores a target by the norm of its projection onto
the range of Y = [Y_1..Y_q] (M x q), using the left singular vectors of Y above a
relative cutoff. The reliability floor is a HEURISTIC: the perturbation dY is measured
against the float64 operand, which the attacker NEVER sees (it is an offline construct
for the measurement). In a deployed attack the floor would be replaced by a
perturbation estimate built from the printed precision itself. The reported rank and
AUROC are therefore descriptive counts under one plausible floor, not a certified
usable-rank theorem.

For each precision (float64, float32, 6 dp, 4 dp) we quantize the RELEASED quantities
(V_t, Lambda_t), form Y_prec, and measure dY = ||Y_prec - Y_f64||_F / ||Y_f64||_F (the
export-rounding perturbation). A direction (relative singular value s_i/s_0 of Y_prec)
is RELIABLE iff s_i/s_0 > dY. We then report:
  - geometric rank (fixed 1e-9 cutoff)
  - reliable rank  (dY cutoff)
  - visible_share_reliable (squared-spectrum of C_H captured by the reliable right
    singular vectors)
  - membership AUROC / TPR@5%FPR on the geometric set vs the reliable set

Output: results/raw/s69_reliable_dirs_chr<C>.json
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
from common import auroc, tpr_at_fpr
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


def _transform(tag):
    if tag == "float64":
        return lambda A: A
    if tag == "float32":
        return lambda A: A.astype(np.float32).astype(np.float64)
    if tag == "round6":
        return lambda A: np.round(A, 6)
    if tag == "round4":
        return lambda A: np.round(A, 4)
    raise ValueError(tag)


def _membership_from_B(B, honest, nonmem, mu, scale, y, cal):
    X = s61.zscore(np.vstack([honest, nonmem]), mu, scale)
    d = np.linalg.norm(X @ B, axis=1) if B.shape[1] else np.zeros(len(X))
    h = np.asarray(cal, bool)
    sgn = 1.0 if auroc(d[h], y[h]) >= 0.5 else -1.0
    return float(auroc(sgn * d[~h], y[~h])), float(tpr_at_fpr(sgn * d[~h], y[~h], 0.05))


def one_split(G, groups, seed, k=s61.K_TOP, block_len=s61.BLOCK_LEN, caps=None):
    roles = s61.split_roles(groups, caps or s61.CAPS, seed)
    honest, partner = G[roles["honest"]], G[roles["partner"]]
    nonmem = G[roles["nonmember"]]
    mu, scale = s61.joint_baseline(honest, partner)
    Zh = s61.zscore(honest, mu, scale)
    C_H = Zh.T @ Zh
    C_Hn = float(np.linalg.norm(C_H))
    n_blocks = G.shape[1] // block_len
    y = np.concatenate([np.ones(len(honest)), np.zeros(len(nonmem))])
    cal = np.random.default_rng(seed + 31).random(len(y)) < 0.5

    Zp0 = s61.zscore(partner, mu, scale)
    lam0, V0 = s61.topk([Zh, Zp0], k)
    K0 = Zp0.T @ Zp0
    rng = np.random.default_rng(seed + 101)
    # released quantities per release, plus the partner's own scatter (known exactly)
    chunks = [(V0, lam0, K0)]
    for t in range(1, max(s61.QS)):
        own, _ = s61.swap_block(partner, rng.choice(partner.shape[0], 2, replace=False),
                                int(rng.integers(0, n_blocks)), block_len)
        Zp = s61.zscore(own, mu, scale)
        lam, V = s61.topk([Zh, Zp], k)
        chunks.append((V, lam, Zp.T @ Zp))
    # the attack operand Y_t = V_t Lambda_t - K_t V_t = C_H V_t (s61.metrics' path):
    # the released V/Λ are quantized, the partner's own K is known exactly.
    Y_f64 = np.concatenate([V * lam[None, :] - K @ V for (V, lam, K) in chunks], axis=1)
    Yf_n = float(np.linalg.norm(Y_f64))

    rows = {}
    for tag in ("float64", "float32", "round6", "round4"):
        T = _transform(tag)
        Ys = [T(V) * T(lam)[None, :] - K @ T(V) for (V, lam, K) in chunks]
        Y = np.concatenate(Ys, axis=1)
        dY = float(np.linalg.norm(Y - Y_f64) / max(1e-12, Yf_n))    # export perturbation
        Wy, sy, Vty = np.linalg.svd(Y, full_matrices=False)
        sy_rel = sy / sy[0]
        floor = max(1e-9, dY)
        keep_geo = sy_rel > 1e-9
        keep_rel = sy_rel > floor
        B_geo = Wy[:, keep_geo]
        B_rel = Wy[:, keep_rel]
        # visible share of C_H captured by the retained left singular directions of Y
        vis_geo = float(np.linalg.norm(C_H @ Wy[:, keep_geo]) / C_Hn) if keep_geo.any() else 0.0
        vis_rel = float(np.linalg.norm(C_H @ Wy[:, keep_rel]) / C_Hn) if keep_rel.any() else 0.0
        auc_geo, tpr_geo = _membership_from_B(B_geo, honest, nonmem, mu, scale, y, cal)
        auc_rel, tpr_rel = _membership_from_B(B_rel, honest, nonmem, mu, scale, y, cal)
        rows[tag] = {"export_perturbation_dY": dY, "reliability_floor": float(floor),
                     "rank_geometric": int(keep_geo.sum()),
                     "rank_reliable": int(keep_rel.sum()),
                     "visible_share_geometric": vis_geo,
                     "visible_share_reliable": vis_rel,
                     "auc_geometric": auc_geo, "tpr05_geometric": tpr_geo,
                     "auc_reliable": auc_rel, "tpr05_reliable": tpr_rel}
    return rows


def run(chrom=22, seeds=5):
    out = {"chrom": chrom, "seeds": seeds, "cohorts": {}}
    for name, (G, groups) in cohorts_for(chrom).items():
        rows = [one_split(G, groups, s) for s in range(seeds)]
        blk = {}
        for tag in ("float64", "float32", "round6", "round4"):
            blk[tag] = {k: float(np.mean([r[tag][k] for r in rows])) for k in rows[0][tag]}
        out["cohorts"][name] = blk
        print("### %s" % name, flush=True)
        for tag in ("float64", "float32", "round6", "round4"):
            b = blk[tag]
            print("   %-8s dY %.2e geom %.1f rel %.1f vis_rel %.3f  auc geo %.3f rel %.3f" % (
                tag, b["export_perturbation_dY"], b["rank_geometric"], b["rank_reliable"],
                b["visible_share_reliable"], b["auc_geometric"], b["auc_reliable"]), flush=True)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s69_reliable_dirs_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
