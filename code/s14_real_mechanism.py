"""S14 - real-data held-out-site mechanism test (S9d analog on 1000G chr22).

Round-2 verdict's cheapest, biggest residual gap: the mechanism claim
(in-sample-ness of the study-fitted basis) is established only in simulation.
Here we reproduce the S9d test on real data: split the real panel into disjoint
sites A and B per super-population plus a shared non-member set, and compare
(a) a site's OWN fitted basis separating its own members vs (b) the OTHER site's
basis separating them (out-of-sample) vs (c) a swapped null.

Output: results/raw/s14_real_mechanism.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import allele_freq, auroc, pca
from real_panel import get_real_panel, SUPERPOP

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def pc_dist(G, loadings, mean):
    X = G.astype(np.float64) - mean[None, :]
    proj = (X @ loadings.T) / np.sqrt(max(1, G.shape[1]))
    return np.linalg.norm(proj, axis=1)


def cal(a):
    return max(float(a), 1 - float(a)), float(a)


def four_way_split(G, sup, seed, per=150):
    """Return A, B, C (two study sites + shared non-member) per super-pop."""
    rng = np.random.default_rng(seed)
    a_i, b_i, c_i = [], [], []
    for p in SUPERPOP:
        idx = rng.permutation(np.where(sup == p)[0])
        per_p = min(per, idx.size // 4)
        a_i.append(idx[:per_p]); b_i.append(idx[per_p:2 * per_p])
        c_i.append(idx[2 * per_p:3 * per_p])
    return (np.concatenate(a_i), np.concatenate(b_i), np.concatenate(c_i))


def run(n_splits=20):
    G, sup, pos = get_real_panel()
    rows = []
    for s in range(n_splits):
        ai, bi, ci = four_way_split(G, sup, s)
        A, B, C = G[ai], G[bi], G[ci]
        _, _, LA = pca(A, k=10); mA = A.mean(0)
        _, _, LB = pca(B, k=10); mB = B.mean(0)
        out = {}
        for name, panel, other, own_L, own_m in (
                ("A", A, B, LA, mA), ("B", B, A, LB, mB)):
            y = np.concatenate([np.ones(len(panel)), np.zeros(len(C))])
            d = pc_dist(panel, own_L, own_m); dn = pc_dist(C, own_L, own_m)
            a = auroc(np.concatenate([d, dn]), y)
            out[f"{name}_in_sample_calibrated"], out[f"{name}_in_sample_naive"] = cal(a)
            o = pc_dist(panel, (LB, LA)[name == "B"], (mB, mA)[name == "B"])
            on = pc_dist(C, (LB, LA)[name == "B"], (mB, mA)[name == "B"])
            a2 = auroc(np.concatenate([o, on]), y)
            out[f"{name}_out_of_sample_calibrated"], out[f"{name}_out_of_sample_naive"] = cal(a2)
            x = pc_dist(other, own_L, own_m); xn = pc_dist(C, own_L, own_m)
            a3 = auroc(np.concatenate([x, xn]), y)
            out[f"{name}_swapped_null_calibrated"] = cal(a3)[0]
        out["seed"] = s
        rows.append(out)
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v)
                          for k, v in out.items()}), flush=True)
    keys = ["A_in_sample_calibrated", "A_out_of_sample_calibrated",
            "A_swapped_null_calibrated", "B_in_sample_calibrated",
            "B_out_of_sample_calibrated", "B_swapped_null_calibrated"]
    agg = {}
    for k in keys:
        agg[k] = float(np.mean([r[k] for r in rows]))
        agg[k + "_std"] = float(np.std([r[k] for r in rows]))
    res = dict(per_seed=rows, agg=agg, n_splits=n_splits,
               meta=dict(n_snps=int(G.shape[1]), n_splits=n_splits))
    with open(os.path.join(OUT, "s14_real_mechanism.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("=== agg ===", json.dumps({k: round(v, 4) for k, v in agg.items()}), flush=True)
    return res


if __name__ == "__main__":
    run()
