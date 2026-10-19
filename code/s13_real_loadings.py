"""S13 - PC-spectrum granularity on REAL 1000G chr22.

The simulation S9c says the released PC leak needs the full loading matrix:
eigenvalue-only release is at chance (0.51) while full loadings give 0.88.
The result-to-claim verdict's central open question is whether the PC channel
exists on real data at all (naive 0.347 / flipped 0.65, S9). This experiment
runs the granularity map on the real panel, so "is the PC channel real?"
becomes "real release of WHICH granularity, under WHICH orientation, at WHAT
magnitude" -- the question governance actually has to answer.

Output: results/raw/s13_real_loadings.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import allele_freq, af_lrt_attack, auroc, pca
from real_panel import get_real_panel, three_way_split

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def pc_dist(G, loadings, mean, k=None):
    L = loadings if k is None else loadings[:k]
    X = G.astype(np.float64) - mean[None, :]
    proj = (X @ L.T) / np.sqrt(max(1, G.shape[1]))
    return np.linalg.norm(proj, axis=1)


def one_split(G, sup, seed):
    sidx, ridx, nidx = three_way_split(sup, per=150, seed=seed)
    site, ref, non = G[sidx], G[ridx], G[nidx]
    _, _, load = pca(site, k=10)
    mean = site.mean(0)
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    out = {}

    def cal(k=None, basis=None):
        L = load if basis is None else basis
        dm = pc_dist(site, L, mean, k=k)
        dn = pc_dist(non, L, mean, k=k)
        a = float(auroc(np.concatenate([-dm, -dn]), y))
        return a, max(a, 1 - a)

    out["full_loadings"], out["full_loadings_calibrated"] = cal()
    for K in (1, 2, 3, 5):
        out[f"top{K}_loadings"], out[f"top{K}_loadings_calibrated"] = cal(k=K)
    rng = np.random.default_rng(seed)
    Lrand = rng.normal(size=load.shape)
    Lrand, _ = np.linalg.qr(Lrand.T); Lrand = Lrand.T
    out["eigenvalues_only"], out["eigenvalues_only_calibrated"] = cal(basis=Lrand)
    out["af_only"] = float(auroc(np.concatenate([
        af_lrt_attack(site, allele_freq(site), allele_freq(ref)),
        af_lrt_attack(non, allele_freq(site), allele_freq(ref))]), y))
    out["seed"] = seed
    return out


def run(split_seeds=tuple(range(20))):
    G, sup, pos = get_real_panel()
    rows = [one_split(G, sup, s) for s in split_seeds]
    for r in rows:
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v)
                          for k, v in r.items()}), flush=True)
    keys = ["full_loadings_calibrated", "top1_loadings_calibrated",
            "top2_loadings_calibrated", "top3_loadings_calibrated",
            "top5_loadings_calibrated", "eigenvalues_only_calibrated",
            "full_loadings", "af_only"]
    agg = {}
    for k in keys:
        agg[k] = float(np.mean([r[k] for r in rows]))
        agg[k + "_std"] = float(np.std([r[k] for r in rows]))
    res = dict(per_seed=rows, agg=agg, n_snps=int(G.shape[1]), n_splits=len(split_seeds))
    with open(os.path.join(OUT, "s13_real_loadings.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("=== agg ===", json.dumps({k: round(v, 4) for k, v in agg.items()}), flush=True)
    return res


if __name__ == "__main__":
    run()
