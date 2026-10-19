"""S35 - AF-dominance frontier shape (companion to the s27 boundary grid).

The s27 grid shows the PC channel loses to the AF channel throughout the
biobank-scale region, but the AF channel there sits at AUROC ~1.0 (saturated).
This block asks the sharper question directly: at what AF-AUROC operating point
does the PC channel stop adding anything? We sweep marker count densely at small
study size (where PC still competes) and locate the crossing, reporting
(af, pc, delta) vs M together with the top-10 variance fraction (a scale-free
proxy for how much between-group structure the fitted spectrum captures).

Output: results/raw/s35_frontier.json
"""
import json
import os
import sys
import time
import numpy as np
import torch

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
import common_gpu as cg
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

OUT = "/data/lvvyongji1/GWAS/results/raw"
DEV = torch.device("cuda")


def cell(nsp, M, fst, seed, k=10, blk=50000):
    pop = cg.simulate_population_host(3, nsp * 3, M, fst, 100, seed=seed, pct="gauss")
    site, ref, non = cg.split_population3_host(pop, nsp, nsp, nsp, seed=seed)
    n_s = site["G"].shape[0]; n_n = non["G"].shape[0]
    ev, L = cg.gram_svd_host(site["G"], k=k, block=blk)
    mu = site["G"].mean(0, dtype=np.float64)
    af_s = site["G"].mean(0, dtype=np.float64) / 2.0
    af_r = ref["G"].mean(0, dtype=np.float64) / 2.0
    am = cg.af_lrt_host(site["G"], af_s, af_r, block=blk)
    an = cg.af_lrt_host(non["G"], af_s, af_r, block=blk)
    dm = cg.proj_dist_host(site["G"], L, mu, block=blk)
    dn = cg.proj_dist_host(non["G"], L, mu, block=blk)
    rng = np.random.default_rng(seed + 313)
    tsm = torch.as_tensor(rng.random(n_s) < 0.5, device=DEV)
    tsn = torch.as_tensor(rng.random(n_n) < 0.5, device=DEV)
    ytr = torch.cat([torch.ones(int(tsm.sum()), device=DEV), torch.zeros(int(tsn.sum()), device=DEV)])
    yte = torch.cat([torch.ones(int((~tsm).sum()), device=DEV), torch.zeros(int((~tsn).sum()), device=DEV)])
    pc_tr = torch.cat([dm[tsm], dn[tsn]])
    sgn = 1.0 if cg.auroc_gpu(-pc_tr, ytr) >= cg.auroc_gpu(pc_tr, ytr) else -1.0
    af_te = torch.cat([am[~tsm], an[~tsn]]); pc_te = sgn * (-torch.cat([dm[~tsm], dn[~tsn]]))
    af_auc = float(cg.auroc_gpu(af_te, yte)); pc_auc = float(cg.auroc_gpu(pc_te, yte))
    ytr_np = ytr.cpu().numpy(); yte_np = yte.cpu().numpy()
    Xtr = np.column_stack([torch.cat([am[tsm], an[tsn]]).cpu().numpy(),
                           torch.cat([dm[tsm], dn[tsn]]).cpu().numpy() * (-sgn)])
    Xte = np.column_stack([af_te.cpu().numpy(), pc_te.cpu().numpy()])
    comb = float(roc_auc_score(yte_np, LogisticRegression(max_iter=1000).fit(Xtr, ytr_np)
                               .predict_proba(Xte)[:, 1]))
    res = dict(af=af_auc, pc=pc_auc, combined=comb, delta=pc_auc - af_auc,
               af_sat=float(af_auc > 0.99), top10_varfrac=float(ev.sum() / (n_s * M)))
    del L, site, ref, non, pop
    torch.cuda.empty_cache()
    return res


def run():
    cells = [(500, M, 0.08) for M in (5000, 10000, 20000, 50000, 100000, 300000, 1000000)]
    cells += [(2000, M, 0.08) for M in (5000, 10000, 20000, 50000, 100000)]
    grid = {}
    t0 = time.time()
    for (nsp, M, fst) in cells:
        rows = []
        for s in range(6):
            try:
                rows.append(cell(nsp, M, fst, s))
            except RuntimeError as e:
                torch.cuda.empty_cache()
                if "out of memory" in str(e).lower():
                    break
                raise
        if not rows:
            continue
        key = f"nsp{nsp}_m{M}"
        agg = {kk: float(np.mean([r[kk] for r in rows])) for kk in rows[0]}
        agg["n_seeds"] = len(rows)
        grid[key] = agg
        print(key, "af=%.3f pc=%.3f delta=%+.3f comb=%.3f" % (
            agg["af"], agg["pc"], agg["delta"], agg["combined"]), flush=True)
        with open(os.path.join(OUT, "s35_frontier.json"), "w") as f:
            json.dump(dict(grid=grid, elapsed_s=round(time.time() - t0, 1)), f, indent=2)
    print("done", round(time.time() - t0, 1), "s", flush=True)


if __name__ == "__main__":
    run()
