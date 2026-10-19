"""S16 - genome-scale real-data replication on independent chromosomes.

The chr22-only real result (PC channel present but BELOW the AF release, at
0.652 vs 0.688, with the in-sample mechanism and the loadings-granularity
ladder both reproducing) is the single-dataset limitation. This runs the same
four studies on chr1 and chr2 (independent chromosomes, ~20x more markers than
chr22) and aggregates, to test whether the real-data ordering replicates.

Studies per chromosome (all disjoint three-way splits, 20 split seeds):
  - frontier: AF tiers (full/2dp/1dp/top200/per-pop) - the S7 ordering
  - orientation: AF vs naive-PC vs calibrated-PC vs AF+PC
  - mechanism: held-out-site (own basis vs other-site basis vs swapped null)
  - granularity: top-k loadings vs eigenvalue-only (calibrated)

Output: results/raw/s16_multichrom.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import allele_freq, af_lrt_attack, auroc, pca
from real_multi import get_real_panel_chr, three_way_split, SUPERPOP

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def cal(a):
    return max(float(a), 1 - float(a))


def pc_dist(G, loadings, mean, k=None):
    L = loadings if k is None else loadings[:k]
    X = G.astype(np.float64) - mean[None, :]
    proj = (X @ L.T) / np.sqrt(max(1, G.shape[1]))
    return np.linalg.norm(proj, axis=1)


def orientation(G, sup, s):
    sidx, ridx, nidx = three_way_split(sup, per=150, seed=s)
    site, ref, non = G[sidx], G[ridx], G[nidx]
    _, _, L = pca(site, k=10); mean = site.mean(0)
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    dm, dn = pc_dist(site, L, mean), pc_dist(non, L, mean)
    af_s, af_r = allele_freq(site), allele_freq(ref)
    am = af_lrt_attack(site, af_s, af_r); an = af_lrt_attack(non, af_s, af_r)
    return dict(
        af_only=float(auroc(np.concatenate([am, an]), y)),
        pc_naive=float(auroc(np.concatenate([dm, dn]), y)),
        pc_calibrated=cal(auroc(np.concatenate([dm, dn]), y)),
        af_plus_pc=None, s=s)


def orientation_combined(G, sup, s):
    """AF+PC combined via logistic regression fit on a train half."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    sidx, ridx, nidx = three_way_split(sup, per=150, seed=s)
    site, ref, non = G[sidx], G[ridx], G[nidx]
    _, _, L = pca(site, k=10); mean = site.mean(0)
    dm, dn = pc_dist(site, L, mean), pc_dist(non, L, mean)
    af_s, af_r = allele_freq(site), allele_freq(ref)
    am = af_lrt_attack(site, af_s, af_r); an = af_lrt_attack(non, af_s, af_r)
    rng = np.random.default_rng(s + 31)
    tsm = rng.random(len(dm)) < 0.5; tsn = rng.random(len(dn)) < 0.5
    ytr = np.concatenate([np.ones(tsm.sum()), np.zeros(tsn.sum())])
    yte = np.concatenate([np.ones((~tsm).sum()), np.zeros((~tsn).sum())])
    pc_tr = np.concatenate([dm[tsm], dn[tsn]])
    sgn = 1.0 if auroc(pc_tr, ytr) >= auroc(-pc_tr, ytr) else -1.0
    Xtr = np.column_stack([np.concatenate([am[tsm], an[tsn]]), sgn * pc_tr])
    Xte = np.column_stack([np.concatenate([am[~tsm], an[~tsn]]),
                           sgn * np.concatenate([dm[~tsm], dn[~tsn]])])
    sc = StandardScaler().fit(Xtr)
    clf = LogisticRegression(max_iter=2000).fit(sc.transform(Xtr), ytr)
    comb = clf.predict_proba(sc.transform(Xte))[:, 1]
    af_te = float(auroc(Xte[:, 0], yte))
    return dict(af_only=af_te, af_plus_pc=float(auroc(comb, yte)),
                pc_calibrated=float(auroc(sgn * Xte[:, 1], yte)), s=s)


def mechanism(G, sup, s):
    rng = np.random.default_rng(s)
    ai, bi, ci = [], [], []
    for p in SUPERPOP:
        idx = rng.permutation(np.where(sup == p)[0])
        per_p = min(150, idx.size // 4)
        ai.append(idx[:per_p]); bi.append(idx[per_p:2 * per_p]); ci.append(idx[2 * per_p:3 * per_p])
    A, B, C = G[np.concatenate(ai)], G[np.concatenate(bi)], G[np.concatenate(ci)]
    _, _, LA = pca(A, k=10); mA = A.mean(0)
    _, _, LB = pca(B, k=10); mB = B.mean(0)
    out = {"s": s}
    for name, panel, oL, oM, othL, othM, other in (
            ("A", A, LA, mA, LB, mB, B), ("B", B, LB, mB, LA, mA, A)):
        y = np.concatenate([np.ones(len(panel)), np.zeros(len(C))])
        d = pc_dist(panel, oL, oM); dn = pc_dist(C, oL, oM)
        out[f"{name}_in"] = cal(auroc(np.concatenate([d, dn]), y))
        oo = pc_dist(panel, othL, othM); on = pc_dist(C, othL, othM)
        out[f"{name}_oos"] = cal(auroc(np.concatenate([oo, on]), y))
        x = pc_dist(other, oL, oM); xn = pc_dist(C, oL, oM)
        out[f"{name}_swapped"] = cal(auroc(np.concatenate([x, xn]), y))
    return out


def granularity(G, sup, s):
    sidx, ridx, nidx = three_way_split(sup, per=150, seed=s)
    site, non = G[sidx], G[nidx]
    _, _, L = pca(site, k=10); mean = site.mean(0)
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    out = {"s": s}
    for K in (1, 2, 3, 5):
        dm = pc_dist(site, L, mean, k=K); dn = pc_dist(non, L, mean, k=K)
        out[f"top{K}"] = cal(auroc(np.concatenate([dm, dn]), y))
    dm = pc_dist(site, L, mean); dn = pc_dist(non, L, mean)
    out["full"] = cal(auroc(np.concatenate([dm, dn]), y))
    rng = np.random.default_rng(s)
    Lr = rng.normal(size=L.shape); Lr, _ = np.linalg.qr(Lr.T); Lr = Lr.T
    dm = pc_dist(site, Lr, mean); dn = pc_dist(non, Lr, mean)
    out["eig_only"] = cal(auroc(np.concatenate([dm, dn]), y))
    return out


def agg(rows, keys):
    return {k: (float(np.mean([r[k] for r in rows])),
                float(np.std([r[k] for r in rows]))) for k in keys}


def run(chroms=(1, 2), n_splits=20):
    res = {"chroms": list(chroms), "n_splits": n_splits, "per_chrom": {}}
    for c in chroms:
        G, sup, pos = get_real_panel_chr(c)
        print(f"=== chr{c} {G.shape} ===", flush=True)
        ori = [orientation_combined(G, sup, s) for s in range(n_splits)]
        mec = [mechanism(G, sup, s) for s in range(n_splits)]
        grn = [granularity(G, sup, s) for s in range(n_splits)]
        res["per_chrom"][f"chr{c}"] = dict(
            n_snps=int(G.shape[1]),
            orientation=agg(ori, ["af_only", "pc_calibrated", "af_plus_pc"]),
            mechanism=agg(mec, ["A_in", "A_oos", "A_swapped", "B_in", "B_oos", "B_swapped"]),
            granularity=agg(grn, ["top1", "top2", "top3", "top5", "full", "eig_only"]),
        )
        print(f"chr{c} orientation",
              json.dumps({k: round(v[0], 4) for k, v in res["per_chrom"][f"chr{c}"]["orientation"].items()}), flush=True)
        print(f"chr{c} mechanism",
              json.dumps({k: round(v[0], 4) for k, v in res["per_chrom"][f"chr{c}"]["mechanism"].items()}), flush=True)
        print(f"chr{c} granularity",
              json.dumps({k: round(v[0], 4) for k, v in res["per_chrom"][f"chr{c}"]["granularity"].items()}), flush=True)
    with open(os.path.join(OUT, "s16_multichrom.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
