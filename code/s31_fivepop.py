"""S31 - All five 1000G super-populations (B4, second half).

The paper's real arm uses 3 super-population pseudo-sites (AFR/EAS/EUR). This
closes the "only 3 super-pops" caveat: repeat the orientation + per-site release
attacks using all five 1000G super-populations (adds SAS and AMR), on chr22 and
chr1+chr2. Reports (a) the PC-vs-AF ordering under 5-pop, (b) a paired 5-pop vs
3-pop comparison, (c) per-superpop AF/PC AUROC.

Output: results/raw/s31_fivepop.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from common import allele_freq, af_lrt_attack, auroc, pca
from real_multi import SUPERPOP5 as SUPER_POP5

OUT = "/data/lvvyongji1/GWAS/results/raw"
SUPS5 = ["AFR", "AMR", "EAS", "EUR", "SAS"]
SUPS3 = ["AFR", "EAS", "EUR"]


def cal(a):
    return max(float(a), 1.0 - float(a))


def pc_dist(G, loadings, mean, k=None):
    L = loadings if k is None else loadings[:k]
    X = G.astype(np.float64) - mean[None, :]
    return np.linalg.norm((X @ L.T) / np.sqrt(max(1, G.shape[1])), axis=1)


def orient_and_score(site, ref, non):
    """AF-only + PC (train-half calibrated) AUROC for one site/ref/non triple."""
    _, _, L = pca(site, k=10); mean = site.mean(0)
    af_s, af_r = allele_freq(site), allele_freq(ref)
    am = af_lrt_attack(site, af_s, af_r); an = af_lrt_attack(non, af_s, af_r)
    dm = pc_dist(site, L, mean); dn = pc_dist(non, L, mean)
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    af_auc = cal(auroc(np.concatenate([am, an]), y))
    ypc = cal(auroc(np.concatenate([dm, dn]), y))
    return dict(af=af_auc, pc=ypc, delta=ypc - af_auc)


def split_by_sups(sup, sups, per=150, seed=0):
    """Per-superpop 3-way disjoint split (study/ref/non), restricted to `sups`."""
    rng = np.random.default_rng(seed)
    si, ri, ni = [], [], []
    for p in sups:
        idx = rng.permutation(np.where(sup == p)[0])
        k = min(per, idx.size // 3)
        si.append(idx[:k]); ri.append(idx[k:2 * k]); ni.append(idx[2 * k:3 * k])
    return np.concatenate(si), np.concatenate(ri), np.concatenate(ni)


def run_one(G, sup, sups, n_seeds=6, per=150):
    sups = list(sups)
    sites = []
    for s in range(n_seeds):
        si, ri, ni = split_by_sups(sup, sups, per=per, seed=s)
        r = orient_and_score(G[si], G[ri], G[ni])
        r["n_site"] = int(si.size)
        r["per_sup"] = {}
        # per-superpop AUROC: score within a single superpop's site/ref/non split
        for p in sups:
            sp_si = si[sup[si] == p]; sp_ni = ni[sup[ni] == p]
            sp_ri = ri[sup[ri] == p]
            if sp_si.size < 20 or sp_ni.size < 20:
                continue
            rr = orient_and_score(G[sp_si], G[sp_ri], G[sp_ni])
            r["per_sup"][p] = rr
        sites.append(r)
    af = np.mean([r["af"] for r in sites]); pc = np.mean([r["pc"] for r in sites])
    d = np.array([r["delta"] for r in sites])
    out = dict(af=float(af), pc=float(pc), delta=float(d.mean()),
               delta_std=float(d.std()), n_site=int(np.mean([r["n_site"] for r in sites])),
               per_sup={p: dict(af=float(np.mean([r["per_sup"][p]["af"] for r in sites if p in r["per_sup"]])),
                                pc=float(np.mean([r["per_sup"][p]["pc"] for r in sites if p in r["per_sup"]])))
                        for p in sups if all(p in r["per_sup"] for r in sites)})
    return out


def main():
    from real_multi import get_multi_chrom, get_real_panel_chr
    res = {}
    PER = 100   # matched per-superpop size for the 5-pop vs 3-pop comparison
    # 1. chr22 3-pop (cached)
    try:
        from real_panel import get_real_panel
        G22, sup22, pos22 = get_real_panel()
        res["chr22_3pop"] = run_one(G22, sup22, SUPS3, per=150)
        print("chr22 3pop af=%.3f pc=%.3f n=%d" % (
            res["chr22_3pop"]["af"], res["chr22_3pop"]["pc"], res["chr22_3pop"]["n_site"]), flush=True)
    except Exception as e:
        res["chr22_3pop"] = dict(error=str(e))
    # 2. chr22 5-pop (small parse)
    try:
        G22_5, sup22_5, _ = get_real_panel_chr(22, sups=SUPER_POP5)
        res["chr22_5pop"] = run_one(G22_5, sup22_5, SUPS5, per=PER)
        print("chr22 5pop af=%.3f pc=%.3f n=%d" % (
            res["chr22_5pop"]["af"], res["chr22_5pop"]["pc"], res["chr22_5pop"]["n_site"]), flush=True)
    except Exception as e:
        res["chr22_5pop"] = dict(error=str(e))
    # 3. chr1+chr2 3-pop (cached) and 5-pop (large parse)
    try:
        G12_3, sup12_3, _ = get_multi_chrom((1, 2), sups=SUPS3)
        res["chr12_3pop"] = run_one(G12_3, sup12_3, SUPS3, per=PER)
        print("chr12 3pop af=%.3f pc=%.3f n=%d" % (
            res["chr12_3pop"]["af"], res["chr12_3pop"]["pc"], res["chr12_3pop"]["n_site"]), flush=True)
    except Exception as e:
        res["chr12_3pop"] = dict(error=str(e))
    try:
        G12_5, sup12_5, _ = get_multi_chrom((1, 2), sups=SUPER_POP5)
        res["chr12_5pop"] = run_one(G12_5, sup12_5, SUPS5, per=PER)
        print("chr12 5pop af=%.3f pc=%.3f n=%d" % (
            res["chr12_5pop"]["af"], res["chr12_5pop"]["pc"], res["chr12_5pop"]["n_site"]), flush=True)
    except Exception as e:
        res["chr12_5pop"] = dict(error=str(e))
    res["sups5"] = SUPS5
    with open(os.path.join(OUT, "s31_fivepop.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("done", flush=True)


if __name__ == "__main__":
    main()
