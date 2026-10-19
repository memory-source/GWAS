"""S27 - GPU biobank-scale boundary grid (host-resident panels).

Locates the PC-vs-AF *dominance boundary* in the (study size, marker count,
Fst) space and tests whether it is a small-study phenomenon. Panels are
HOST-resident int8 (up to ~90 GB); the exact PCA Gram and attack scores STREAM
marker blocks to the GPU, so GPU memory is O(n^2), independent of M.

Design note: the PCA/attack cost is set by the STUDY size n = 3*nsp (nsp =
study individuals per group), NOT the pool size -- so we sweep nsp directly.
The pool is nsp*4 per group so the 3-way disjoint split (site/ref/non) has
exactly nsp each.

Per cell: AF-only, train-half-calibrated PC, logistic-combined, paired
PC-AF gap, exact top-10 eigenvalue spectrum with eigenvalue-only and
reference-basis controls (in-sample mechanism), and fst_like.

Usage:  python3 code/s27_gpu_scale.py [shard] [nshards]
Output: results/raw/s27_gpu_scale_shard<s>.json
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


def fst_like_np(G, anc):
    af = G.mean(0, dtype=np.float64) / 2.0
    hs = np.mean([2 * (G[anc == k].mean(0, dtype=np.float64) / 2.0)
                  * (1 - G[anc == k].mean(0, dtype=np.float64) / 2.0)
                  for k in np.unique(anc)], axis=0)
    ht = 2 * af * (1 - af)
    return float(np.nanmean((ht - hs) / np.where(ht > 0, ht, np.nan)))


def auroc_lr(tr_scores, tr_y, te_scores, te_y):
    """Logistic AUROC on the columns of `tr_scores` (shape (n, k)); with k=2 it
    is the joint AF+PC combination, matching s12's combined arm."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    Xtr = np.asarray(tr_scores, float).reshape(len(tr_y), -1)
    Xte = np.asarray(te_scores, float).reshape(len(te_y), -1)
    lr = LogisticRegression(max_iter=1000).fit(Xtr, tr_y)
    return float(roc_auc_score(te_y, lr.predict_proba(Xte)[:, 1]))


def cell(nsp, M, fst, seed, k=10, blk=50000):
    pop = cg.simulate_population_host(n_groups=3, n_per_group=nsp * 3, n_snps=M,
                                      fst=fst, ld_blocks=100, seed=seed, pct="gauss")
    site, ref, non = cg.split_population3_host(pop, nsp, nsp, nsp, seed=seed)
    n_s = site["G"].shape[0]
    ev, L = cg.gram_svd_host(site["G"], k=k, block=blk)      # exact, top-k
    mu = site["G"].mean(0, dtype=np.float64)
    af_s = site["G"].mean(0, dtype=np.float64) / 2.0
    af_r = ref["G"].mean(0, dtype=np.float64) / 2.0

    # ---- AF-LRT arm
    am = cg.af_lrt_host(site["G"], af_s, af_r, block=blk)
    an = cg.af_lrt_host(non["G"], af_s, af_r, block=blk)
    # ---- PC arm: projection distance to study centroid, train-half calibrated
    dm = cg.proj_dist_host(site["G"], L, mu, block=blk)
    dn = cg.proj_dist_host(non["G"], L, mu, block=blk)
    n_n = non["G"].shape[0]

    rng = np.random.default_rng(seed + 313)
    tsm = torch.as_tensor(rng.random(n_s) < 0.5, device=DEV)
    tsn = torch.as_tensor(rng.random(n_n) < 0.5, device=DEV)
    ytr = torch.cat([torch.ones(int(tsm.sum()), device=DEV),
                     torch.zeros(int(tsn.sum()), device=DEV)])
    yte = torch.cat([torch.ones(int((~tsm).sum()), device=DEV),
                     torch.zeros(int((~tsn).sum()), device=DEV)])
    # orient distance so members score HIGHER on the calibrated score
    pc_tr = torch.cat([dm[tsm], dn[tsn]])
    sgn = 1.0 if cg.auroc_gpu(-pc_tr, ytr) >= cg.auroc_gpu(pc_tr, ytr) else -1.0
    pcm = sgn * (-dm); pcn = sgn * (-dn)                     # score: high = member
    af_tr = torch.cat([am[tsm], an[tsn]]); af_te = torch.cat([am[~tsm], an[~tsn]])
    pc_te = torch.cat([pcm[~tsm], pcn[~tsn]])
    af_auc = float(cg.auroc_gpu(af_te, yte))
    pc_auc = float(cg.auroc_gpu(pc_te, yte))
    ytr_np = ytr.cpu().numpy(); yte_np = yte.cpu().numpy()
    pc_tr_te = torch.cat([pcm[tsm], pcn[tsn]])          # train-half PC (oriented)
    comb = auroc_lr(np.column_stack([af_tr.cpu().numpy(), pc_tr_te.cpu().numpy()]),
                    ytr_np,
                    np.column_stack([af_te.cpu().numpy(), pc_te.cpu().numpy()]),
                    yte_np)
    mu_af = float(am.mean()); mu_an = float(an.mean())

    # ---- in-sample mechanism controls: eigenvalue-only (no loadings) and
    #      reference-fitted basis. Both should collapse the PC channel.
    ev_sum = float(ev.sum())
    ev_top10_frac = float(ev.sum() / (n_s * M))              # var-frac scale
    res = dict(af_only=af_auc, pc_calibrated=pc_auc, combined=comb,
               delta_pc_minus_af=pc_auc - af_auc,
               af_member_mean=mu_af, af_non_mean=mu_an,
               ev_top1_frac=float(ev[0] / (n_s * M)),
               ev_top10_frac=ev_top10_frac,
               fst_like=fst_like_np(site["G"], site["ancestry"]),
               n_site=n_s)
    rngb = np.random.default_rng(seed + 999)
    # reference-fitted basis (random half of the reference), centred on the
    # *reference* mean -- the in-sample mechanism control. Should collapse the
    # PC channel toward the structure-only baseline.
    rf = ref["G"]
    nrf = rf.shape[0]
    ridx = rngb.permutation(nrf)[:nrf // 2]
    mu_r = rf[ridx].mean(0, dtype=np.float64)
    _, L_r = cg.gram_svd_host(rf, indices=ridx, k=k, block=blk)
    rdm = cg.proj_dist_host(site["G"], L_r, mu_r, block=blk)
    rdn = cg.proj_dist_host(non["G"], L_r, mu_r, block=blk)
    pc_r_tr = torch.cat([rdm[tsm], rdn[tsn]])
    sgn_r = 1.0 if cg.auroc_gpu(-pc_r_tr, ytr) >= cg.auroc_gpu(pc_r_tr, ytr) else -1.0
    rc_te = torch.cat([sgn_r * (-rdm[~tsm]), sgn_r * (-rdn[~tsn])])
    res["pc_refbasis"] = float(cg.auroc_gpu(rc_te, yte))
    res["pc_in_sample_gap"] = pc_auc - res["pc_refbasis"]
    # eigenvalue-only control (random orthonormal basis, same k, no fitted
    # directions): at chance everywhere if the leak is in the loadings. This is
    # a DIFFERENT control from the reference-fitted basis above.
    rng_e = np.random.default_rng(seed + 4242)
    L_e = rng_e.normal(size=(k, M)).astype(np.float32)
    Q, _ = np.linalg.qr(L_e.T)                          # (M, k), orthonormal cols
    L_e = torch.as_tensor(Q.T, device=DEV)              # (k, M)
    edm = cg.proj_dist_host(site["G"], L_e, mu, block=blk)
    edn = cg.proj_dist_host(non["G"], L_e, mu, block=blk)
    pc_e_tr = torch.cat([edm[tsm], edn[tsn]])
    sgn_e = 1.0 if cg.auroc_gpu(-pc_e_tr, ytr) >= cg.auroc_gpu(pc_e_tr, ytr) else -1.0
    ec_te = torch.cat([sgn_e * (-edm[~tsm]), sgn_e * (-edn[~tsn])])
    res["pc_eigonly"] = float(cg.auroc_gpu(ec_te, yte))
    del L, L_r, L_e, site, ref, non, pop
    torch.cuda.empty_cache()
    return res


def grid_config():
    # cost-aware (nsp, M, fst) cells. The simulator draw count is ~3*nsp*M per
    # cell (measured numpy PCG64 ~2e8 draws/s), so large-n cells use smaller M
    # to keep every cell under ~2.5e9 draws. Spans the (study size) x (marker
    # count) corners needed to locate the PC-dominance boundary.
    cells = []
    for fst in (0.04, 0.08):
        for nsp, M in [(500, 100000), (500, 1000000),
                       (2000, 100000), (2000, 300000),
                       (8000, 50000), (8000, 100000)]:
            cells.append((nsp, M, fst))
    return dict(cells=cells, n_seeds=3, budget=2e9)


def run(shard=0, nshards=1):
    cfg = grid_config()
    combos = cfg["cells"][shard::nshards]
    out_path = os.path.join(OUT, f"s27_gpu_scale_shard{shard}.json")
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
        agg = {}
        for kk in rows[0]:
            agg[kk] = float(np.mean([r[kk] for r in rows]))
            if kk not in ("n_site",):
                agg[kk + "_std"] = float(np.std([r[kk] for r in rows]))
        agg["n_seeds_used"] = len(rows)
        grid[key] = agg
        print(key, json.dumps({kk: round(agg[kk], 4) for kk in
                               ("af_only", "pc_calibrated", "combined",
                                "delta_pc_minus_af", "pc_refbasis", "pc_eigonly",
                                "pc_in_sample_gap", "fst_like")
                               if kk in agg}), flush=True)
        with open(out_path, "w") as f:
            json.dump(dict(grid=grid, cfg=cfg, shard=shard, nshards=nshards,
                           elapsed_s=round(time.time() - t0, 1)), f, indent=2)
    print(f"shard{shard} done", round(time.time() - t0, 1), "s", flush=True)


if __name__ == "__main__":
    sh = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    ns = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    run(sh, ns)
