"""S28 - Batched permutation nulls + bootstrap CIs (B3).

The paper's frontier numbers carry 300-draw/5-seed CIs. This block replaces
them with (a) a 1e5 stratified label-permutation null for the AF and PC AUROC
on each key cell -- validating that neither channel's discrimination is a
structure artifact -- and (b) bootstrap CIs on the paired PC-AF gap.

The permutation kernel exploits the rank identity for AUROC: with scores fixed,
AUROC under a permuted label vector = (rank-sum of positives - npos(npos+1)/2)
/ (npos*nneg). Permuting labels is equivalent to permuting the (fixed) ranks,
so 1e5 permutations are one batched gather + reduction. Stratification: labels
are permuted WITHIN each ancestry group (preserving the group-wise member
counts), which is the null "no member/non-member information beyond ancestry".

Usage:  python3 code/s28_permnull.py [shard] [nshards]
Output: results/raw/s28_permnull_shard<s>.json
"""
import json
import os
import sys
import time
import numpy as np
import torch

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
import common_gpu as cg

OUT = "/data/lvvyongji1/GWAS/results/raw"
DEV = torch.device("cuda")
NPERM = 20000


def rank_auroc_batched(scores, label_batch):
    """scores (N,), label_batch (P, N) float -> (P,) AUROC via rank-sum."""
    s = scores.flatten()
    order = torch.argsort(s)
    ranks = torch.empty_like(s)
    ranks[order] = torch.arange(1, s.numel() + 1, device=s.device, dtype=torch.float32)
    npos = label_batch.sum(1)
    nneg = s.numel() - npos
    sr = label_batch @ ranks                       # (P,)
    return (sr - npos * (npos + 1) / 2) / (npos * nneg)


def stratified_perms(y_np, anc_np, nperm, seed, chunk=100000):
    """Yield chunks of (p, N) permuted label rows, permuting WITHIN ancestry groups.

    Vectorized: for each stratum of size n_g with k_g positives, draw p random
    scores and mark the k_g lowest-scoring positions as positive -- a uniform
    draw (without replacement) of which positions are positive."""
    rng = np.random.default_rng(seed)
    groups = [np.where(anc_np == g)[0] for g in np.unique(anc_np)]
    kpos = [int(y_np[idx].sum()) for idx in groups]
    N = len(y_np)
    done = 0
    while done < nperm:
        p = min(chunk, nperm - done)
        lab = np.zeros((p, N), dtype=np.float32)
        for idx, kg in zip(groups, kpos):
            ng = len(idx)
            if kg == 0:
                continue
            if kg == ng:
                lab[:, idx] = 1.0
                continue
            order = np.argsort(rng.random((p, ng)), axis=1)[:, :kg]   # (p, kg)
            lab[np.arange(p)[:, None], idx[order]] = 1.0
        yield torch.as_tensor(lab, device=DEV)
        done += p


def null_pvalue(scores, y, anc, nperm=NPERM, seed=0):
    """P(permuted AUROC >= observed) plus the null mean/std."""
    obs = float(cg.auroc_gpu(scores, torch.as_tensor(y, device=DEV, dtype=torch.float32)))
    ge = 0
    tot = 0
    vals = []
    for lab in stratified_perms(y.astype(np.float32), anc, nperm, seed):
        a = rank_auroc_batched(scores, lab)
        vals.append(a)
        ge += int((a >= obs).sum().item())
        tot += a.numel()
    allv = torch.cat(vals)
    return dict(observed=obs, perm_mean=float(allv.mean()), perm_std=float(allv.std()),
                p_value=(ge + 1) / (tot + 1), n_perm=tot)


def bootstrap_ci(delta_np, n_boot=10000, seed=0):
    """Bootstrap CI over per-seed paired deltas (PC-AF)."""
    rng = np.random.default_rng(seed)
    d = np.asarray(delta_np)
    if d.size == 1:
        return dict(mean=float(d[0]), lo=float(d[0]), hi=float(d[0]), p_gt0=float(d[0] > 0))
    idx = rng.integers(0, d.size, size=(n_boot, d.size))
    boot = d[idx].mean(1)
    return dict(mean=float(d.mean()),
                lo=float(np.percentile(boot, 2.5)),
                hi=float(np.percentile(boot, 97.5)),
                p_gt0=float((d > 0).mean()))


def cell(nsp, M, fst, seed, k=10, blk=50000, nperm=NPERM):
    pop = cg.simulate_population_host(n_groups=3, n_per_group=nsp * 3, n_snps=M,
                                      fst=fst, ld_blocks=100, seed=seed, pct="gauss")
    site, ref, non = cg.split_population3_host(pop, nsp, nsp, nsp, seed=seed)
    ev, L = cg.gram_svd_host(site["G"], k=k, block=blk)
    mu = site["G"].mean(0, dtype=np.float64)
    af_s = site["G"].mean(0, dtype=np.float64) / 2.0
    af_r = ref["G"].mean(0, dtype=np.float64) / 2.0
    am = cg.af_lrt_host(site["G"], af_s, af_r, block=blk)
    an = cg.af_lrt_host(non["G"], af_s, af_r, block=blk)
    dm = cg.proj_dist_host(site["G"], L, mu, block=blk)
    dn = cg.proj_dist_host(non["G"], L, mu, block=blk)
    af = torch.cat([am, an])
    # orient the PC score on a held-out half (distance sign is not fixed a priori:
    # members can be *closer to or farther from* the study centroid depending on
    # the regime), then test the oriented score on the full sample.
    n_s, n_n = site["G"].shape[0], non["G"].shape[0]
    rng = np.random.default_rng(seed + 55)
    tsm = torch.as_tensor(rng.random(n_s) < 0.5, device=DEV)
    tsn = torch.as_tensor(rng.random(n_n) < 0.5, device=DEV)
    ytr = torch.cat([torch.ones(int(tsm.sum()), device=DEV),
                     torch.zeros(int(tsn.sum()), device=DEV)])
    dtr = torch.cat([dm[tsm], dn[tsn]])
    sgn = 1.0 if cg.auroc_gpu(-dtr, ytr) >= cg.auroc_gpu(dtr, ytr) else -1.0
    pc = sgn * (-torch.cat([dm, dn]))                # score: high = member
    y = np.r_[np.ones(n_s), np.zeros(n_n)]
    anc = np.r_[site["ancestry"], non["ancestry"]]
    afn = null_pvalue(af, y, anc, nperm=nperm, seed=seed + 1)
    pcn = null_pvalue(pc, y, anc, nperm=nperm, seed=seed + 2)
    af_obs = cg.auroc_gpu(af, torch.as_tensor(y, device=DEV, dtype=torch.float32))
    pc_obs = cg.auroc_gpu(pc, torch.as_tensor(y, device=DEV, dtype=torch.float32))
    out = dict(af=afn, pc=pcn,
               delta_pc_minus_af=float(pc_obs - af_obs))
    del L, site, ref, non, pop
    torch.cuda.empty_cache()
    return out


def grid_config():
    # small-study cells only for the permutation nulls (they are the expensive
    # part and the claim to nail is the small-n regime). n=8000 cells are left to
    # s27's boundary grid.
    cells = [(500, 100000, 0.08), (500, 500000, 0.08),
             (2000, 50000, 0.08), (500, 100000, 0.04), (2000, 50000, 0.04)]
    return dict(cells=cells, n_seeds=3, nperm=NPERM)


def run(shard=0, nshards=1):
    cfg = grid_config()
    combos = cfg["cells"][shard::nshards]
    out_path = os.path.join(OUT, f"s28_permnull_shard{shard}.json")
    grid = {}
    t0 = time.time()
    for (nsp, M, fst) in combos:
        rows = []
        for s in range(cfg["n_seeds"]):
            try:
                rows.append(cell(nsp, M, fst, s))
            except RuntimeError as e:
                torch.cuda.empty_cache()
                if "out of memory" in str(e).lower():
                    print("OOM", nsp, M, fst, flush=True)
                    break
                raise
            finally:
                import gc
                gc.collect()
        if not rows:
            continue
        key = f"nsp{nsp}_m{M}_f{fst}"
        deltas = [r["delta_pc_minus_af"] for r in rows]
        grid[key] = dict(
            af_auc=float(np.mean([r["af"]["observed"] for r in rows])),
            pc_auc=float(np.mean([r["pc"]["observed"] for r in rows])),
            af_perm_mean=float(np.mean([r["af"]["perm_mean"] for r in rows])),
            pc_perm_mean=float(np.mean([r["pc"]["perm_mean"] for r in rows])),
            af_p_max=float(max(r["af"]["p_value"] for r in rows)),
            pc_p_max=float(max(r["pc"]["p_value"] for r in rows)),
            delta=float(np.mean(deltas)),
            delta_ci=bootstrap_ci(deltas, seed=1),
            n_seeds_used=len(rows), n_perm=cfg["nperm"])
        g = grid[key]
        print(key, "af=%.3f(p~%.0e) pc=%.3f(p~%.0e) delta=%.3f [%.3f,%.3f]" % (
            g["af_auc"], g["af_p_max"], g["pc_auc"], g["pc_p_max"],
            g["delta"], g["delta_ci"]["lo"], g["delta_ci"]["hi"]), flush=True)
        with open(out_path, "w") as f:
            json.dump(dict(grid=grid, cfg=cfg, shard=shard, nshards=nshards,
                           elapsed_s=round(time.time() - t0, 1)), f, indent=2)
    print(f"shard{shard} done", round(time.time() - t0, 1), "s", flush=True)


if __name__ == "__main__":
    sh = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    ns = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    run(sh, ns)
