"""S30 - Real marker-density robustness (B4, first half).

The paper's real arm thins 1000G to ~5k markers/chr via a 15 kb distance window.
Reviewers flag that. This re-runs the orientation attack at denser marker tiers
to test whether the PC-vs-AF ordering is a density artifact. Tiers (per chr):
  D_thin  : 15 kb window (the paper's value)
  D_med   : 3 kb window
  D_dense : 300 bp window (near all common SNPs)
on chr22 and chr1+chr2. Reports AF/PC/combined AUROC vs density.

Output: results/raw/s30_density.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from common import allele_freq, af_lrt_attack, auroc, pca

OUT = "/data/lvvyongji1/GWAS/results/raw"
SUPS3 = ("AFR", "EAS", "EUR")
TIERS = {"D_thin": 15000, "D_med": 3000, "D_dense": 300}


def cal(a):
    return max(float(a), 1.0 - float(a))


def pc_dist(G, loadings, mean, k=None):
    L = loadings if k is None else loadings[:k]
    X = G.astype(np.float64) - mean[None, :]
    return np.linalg.norm((X @ L.T) / np.sqrt(max(1, G.shape[1])), axis=1)


def split_by_sup(sup, per=150, seed=0):
    rng = np.random.default_rng(seed)
    si, ri, ni = [], [], []
    for p in SUPS3:
        idx = rng.permutation(np.where(sup == p)[0])
        k = min(per, idx.size // 3)
        si.append(idx[:k]); ri.append(idx[k:2 * k]); ni.append(idx[2 * k:3 * k])
    return np.concatenate(si), np.concatenate(ri), np.concatenate(ni)


def one_tier(G, sup, n_seeds=6):
    afs, pcs, ds = [], [], []
    for s in range(n_seeds):
        si, ri, ni = split_by_sup(sup, seed=s)
        site, ref, non = G[si], G[ri], G[ni]
        _, _, L = pca(site, k=10); mean = site.mean(0)
        am = af_lrt_attack(site, allele_freq(site), allele_freq(ref))
        an = af_lrt_attack(non, allele_freq(site), allele_freq(ref))
        dm = pc_dist(site, L, mean); dn = pc_dist(non, L, mean)
        y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
        af = cal(auroc(np.concatenate([am, an]), y))
        pc = cal(auroc(np.concatenate([dm, dn]), y))
        afs.append(af); pcs.append(pc); ds.append(pc - af)
    return dict(af=float(np.mean(afs)), pc=float(np.mean(pcs)),
                delta=float(np.mean(ds)), delta_std=float(np.std(ds)),
                n_markers=int(G.shape[1]), n_seeds=n_seeds)


def main():
    from real_multi import get_real_panel_chr, get_multi_chrom
    res = {}
    for c in (22, 1, 2):
        res[f"chr{c}"] = {}
        for name, w in TIERS.items():
            try:
                G, sup, pos = get_real_panel_chr(c, ld_window=w, max_snps=500000)
                r = one_tier(G, sup)
                res[f"chr{c}"][name] = r
                print(f"chr{c} {name} m={r['n_markers']} af={r['af']:.3f} pc={r['pc']:.3f} d={r['delta']:+.3f}",
                      flush=True)
            except Exception as e:
                res[f"chr{c}"][name] = dict(error=str(e))
    with open(os.path.join(OUT, "s30_density.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("done", flush=True)


if __name__ == "__main__":
    main()
