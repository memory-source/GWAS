"""S17 - larger, more realistic simulation regime.

Tests whether the PC-dominance and the in-sample mechanism survive at realistic
site sizes:
  - site size sweep: n_per_group in {1200, 3000} at M=10000 SNPs, ld_blocks=100
  - one zero-LD cell (ld_blocks=0) at n=3000, M=10000
  - 8 seeds per cell
Reported per cell: AF-only, calibrated PC (train-half), combined; paired PC-AF
with std.

Output: results/raw/s17_scale.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def cell(n_per_group, n_snps, ld_blocks, n_seeds):
    afv, pcv, combv = [], [], []
    for seed in range(n_seeds):
        pop = simulate_population(n_groups=3, n_per_group=n_per_group,
                                  n_snps=n_snps, ld_blocks=ld_blocks, seed=seed)
        per = max(50, n_per_group // 3)
        site, ref, non = split_population3(pop, per, per, per)
        rel = release_objects(site)
        mask = rel["mask"]
        if mask.sum() < 10:
            continue
        af_s, af_r = rel["af"], allele_freq(ref["G"])
        mean_s = site["G"][:, mask].mean(0)
        L = rel["loadings"]

        def pc(G):
            X = G[:, mask].astype(np.float64) - mean_s[None, :]
            return np.linalg.norm((X @ L[:, mask].T) / np.sqrt(mask.sum()), axis=1)

        dm, dn = pc(site["G"]), pc(non["G"])
        am = af_lrt_attack(site["G"], af_s, af_r, mask=mask)
        an = af_lrt_attack(non["G"], af_s, af_r, mask=mask)
        rng = np.random.default_rng(seed + 991)
        tsm = rng.random(len(dm)) < 0.5; tsn = rng.random(len(dn)) < 0.5
        ytr = np.concatenate([np.ones(tsm.sum()), np.zeros(tsn.sum())])
        yte = np.concatenate([np.ones((~tsm).sum()), np.zeros((~tsn).sum())])
        pc_tr = np.concatenate([dm[tsm], dn[tsn]])
        sgn = 1.0 if auroc(pc_tr, ytr) >= auroc(-pc_tr, ytr) else -1.0
        af_tr = np.concatenate([am[tsm], an[tsn]])
        af_te = float(auroc(np.concatenate([am[~tsm], an[~tsn]]), yte))
        pc_te = float(auroc(sgn * np.concatenate([dm[~tsm], dn[~tsn]]), yte))
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        Xtr = np.column_stack([af_tr, sgn * pc_tr])
        Xte = np.column_stack([np.concatenate([am[~tsm], an[~tsn]]),
                              sgn * np.concatenate([dm[~tsm], dn[~tsn]])])
        sc = StandardScaler().fit(Xtr)
        comb = LogisticRegression(max_iter=1000).fit(sc.transform(Xtr), ytr).predict_proba(
            sc.transform(Xte))[:, 1]
        afv.append(af_te); pcv.append(pc_te); combv.append(float(auroc(comb, yte)))
        print(f"  seed {seed} done", flush=True)
    d = np.array(pcv) - np.array(afv)
    return dict(af_only=float(np.mean(afv)), pc_calibrated=float(np.mean(pcv)),
                combined=float(np.mean(combv)),
                delta_pc_minus_af=float(d.mean()), delta_std=float(d.std()),
                delta_comb_minus_af=float((np.array(combv) - np.array(afv)).mean()),
                n=len(afv))


def run():
    out = {"cells": {}}
    t0 = time.time()

    def flush():
        with open(os.path.join(OUT, "s17_scale.json"), "w") as f:
            json.dump(out, f, indent=2)

    for npg in (1200, 3000):
        print(f"--- n_per_group={npg}, M=10000 ---", flush=True)
        c = cell(npg, 10000, 100, 8)
        out["cells"][f"n{npg}_m10000_ld100"] = c
        print(json.dumps({k: round(v, 4) if isinstance(v, float) else v for k, v in c.items()}),
              flush=True)
        print("elapsed", round(time.time() - t0, 1), "s", flush=True)
        flush()
    print("--- ld_blocks=0, n=3000, M=10000 ---", flush=True)
    c = cell(3000, 10000, 0, 8)
    out["cells"]["n3000_m10000_ld0"] = c
    print(json.dumps({k: round(v, 4) if isinstance(v, float) else v for k, v in c.items()}),
          flush=True)
    flush()
    return out


if __name__ == "__main__":
    run()
