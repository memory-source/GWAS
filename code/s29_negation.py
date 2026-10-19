"""S29 - Paired negation sweep + exact k-ladder (B3).

Two tightens on the paper's central ordering claim:

(a) PAIRED NEGATION: the paper reports delta = PC - AF with 5-30 seeds. Here we
    run N=600 seeds at the two operating points that carry the claim (the
    small-study corner where PC > AF, and a saturated-AF corner), computing the
    per-seed paired delta and P(delta > 0) plus a bootstrap CI. This converts
    "delta = +0.10 (std 0.02)" into a sign-tested, tightly-CI'd statement.

(b) EXACT k-LADDER: the paper's granularity arm jumps k in {1,2,3,5}. Here we
    compute the PC-projection AUROC for every k in 1..10 (exact top-k PCA), plus
    an eigenvalue-only control (random orthonormal basis, no fitted directions),
    on three corpora: parametric sim, real 1000G chr1/2/22, and msprime. This
    pins whether the leak saturates at a measured k and that it lives in the
    LOADINGS, not the eigenvalues.

Usage:  python3 code/s29_negation.py [shard] [nshards]
Output: results/raw/s29_negation_shard<s>.json
"""
import json
import os
import sys
import time
import numpy as np

# cap BLAS threads: 12 shards x default 96 threads thrash a shared box
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from common import simulate_population, split_population3, release_objects, \
    allele_freq, af_lrt_attack, auroc, pca

OUT = "/data/lvvyongji1/GWAS/results/raw"


def cal(a):
    return max(float(a), 1.0 - float(a))


def pc_dist_k(G, mask, loadings, mean, k=None):
    L = loadings if k is None else loadings[:k]
    X = G[:, mask].astype(np.float64) - mean[None, :]
    return np.linalg.norm((X @ L[:, mask].T) / np.sqrt(max(1, mask.sum())), axis=1)


# ---------------------------------------------------------------- (a) paired ---
def paired_seed(npg, M, fst, seed, per_frac=0.3):
    pop = simulate_population(n_groups=3, n_per_group=npg, n_snps=M, fst=fst, seed=seed)
    per = max(10, int(npg * per_frac))
    site, ref, non = split_population3(pop, per, per, per)
    rel = release_objects(site)
    mask = rel["mask"]
    if mask.sum() < 5:
        return None
    af_s, af_r = rel["af"], allele_freq(ref["G"])
    mean_s = site["G"][:, mask].mean(0)
    load_s = rel["loadings"]                      # (k, M)
    dm = pc_dist_k(site["G"], mask, load_s, mean_s)
    dn = pc_dist_k(non["G"], mask, load_s, mean_s)
    am = af_lrt_attack(site["G"], af_s, af_r, mask=mask)
    an = af_lrt_attack(non["G"], af_s, af_r, mask=mask)
    rng = np.random.default_rng(seed + 7717)
    tsm = rng.random(len(dm)) < 0.5
    tsn = rng.random(len(dn)) < 0.5
    ytr = np.concatenate([np.ones(tsm.sum()), np.zeros(tsn.sum())])
    yte = np.concatenate([np.ones((~tsm).sum()), np.zeros((~tsn).sum())])
    pc_tr = np.concatenate([dm[tsm], dn[tsn]])
    sgn = 1.0 if auroc(pc_tr, ytr) >= auroc(-pc_tr, ytr) else -1.0
    af_te = np.concatenate([am[~tsm], an[~tsn]])
    pc_te = sgn * np.concatenate([dm[~tsm], dn[~tsn]])
    a_af = auroc(af_te, yte); a_pc = auroc(pc_te, yte)
    return dict(af=float(a_af), pc=float(a_pc), delta=float(a_pc - a_af))


def paired_sweep(cells, n_seeds, seed0=0):
    out = {}
    for (npg, M, fst) in cells:
        rows = []
        for s in range(seed0, seed0 + n_seeds):
            r = paired_seed(npg, M, fst, s)
            if r is not None:
                rows.append(r)
        if not rows:
            continue
        d = np.array([r["delta"] for r in rows])
        af = np.array([r["af"] for r in rows]); pc = np.array([r["pc"] for r in rows])
        rng = np.random.default_rng(7)
        boot = d[rng.integers(0, d.size, size=(10000, d.size))].mean(1)
        key = f"n{npg}_m{M}_f{fst}"
        out[key] = dict(
            n_seeds=len(rows), seed0=seed0,
            af_mean=float(af.mean()), pc_mean=float(pc.mean()),
            delta_mean=float(d.mean()), delta_std=float(d.std()),
            delta_ci_lo=float(np.percentile(boot, 2.5)),
            delta_ci_hi=float(np.percentile(boot, 97.5)),
            p_delta_gt0=float((d > 0).mean()),
            delta_min=float(d.min()), delta_max=float(d.max()),
            # per-seed deltas so shards can be pooled honestly downstream
            per_seed_delta=[float(x) for x in d])
        g = out[key]
        print(key, "delta=%.4f [%.4f,%.4f] P(>0)=%.3f (n=%d)" % (
            g["delta_mean"], g["delta_ci_lo"], g["delta_ci_hi"],
            g["p_delta_gt0"], g["n_seeds"]), flush=True)
    return out


# ---------------------------------------------------------------- (b) ladder ---
def ladder_from_panels(site, ref, non, ks=range(1, 11)):
    """Exact top-10 PCA; PC AUROC at every k; eigenvalue-only control."""
    mask = np.ones(site.shape[1], bool)
    _, _, L = pca(site, k=10)
    mean = site.mean(0)
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    am = af_lrt_attack(site, allele_freq(site), allele_freq(ref))
    an = af_lrt_attack(non, allele_freq(site), allele_freq(ref))
    af_auc = cal(auroc(np.concatenate([am, an]), y))
    lad = {}
    for k in ks:
        dm = pc_dist_k(site, mask, L, mean, k=k)
        dn = pc_dist_k(non, mask, L, mean, k=k)
        lad[k] = cal(auroc(np.concatenate([dm, dn]), y))
    # eigenvalue-only control: a random orthonormal basis (no fitted
    # directions), same k, centred on the study mean.
    rng = np.random.default_rng(0)
    Lr = rng.normal(size=L.shape)
    Lr, _ = np.linalg.qr(Lr.T)
    Lr = Lr.T
    dm = pc_dist_k(site, mask, Lr, mean); dn = pc_dist_k(non, mask, Lr, mean)
    eig_only = cal(auroc(np.concatenate([dm, dn]), y))
    return dict(af_only=af_auc, ladder=lad, eig_only=eig_only,
                n_site=int(len(site)), n_markers=int(site.shape[1]))


def ladder_sim(npg=400, M=2000, seeds=(0, 1, 2, 3, 4)):
    acc = {}
    for s in seeds:
        pop = simulate_population(n_groups=3, n_per_group=npg, n_snps=M, seed=s)
        per = max(10, int(npg * 0.3))
        site, ref, non = split_population3(pop, per, per, per)
        r = ladder_from_panels(site["G"], ref["G"], non["G"])
        acc.setdefault("af_only", []).append(r["af_only"])
        acc.setdefault("eig_only", []).append(r["eig_only"])
        for k, v in r["ladder"].items():
            acc.setdefault(k, []).append(v)
    return {str(k): dict(mean=float(np.mean(v)), std=float(np.std(v)))
            for k, v in acc.items()}


def ladder_real():
    from real_panel import get_real_panel
    from real_multi import get_real_panel_chr
    out = {}
    try:
        G, sup, pos = get_real_panel()
        site, ref, non = _split_by_sup(G, sup)
        out["chr22"] = ladder_from_panels(site, ref, non)
    except Exception as e:
        out["chr22"] = dict(error=str(e))
    for chrom in (1, 2):
        try:
            G, sup, pos = get_real_panel_chr(chrom)
            site, ref, non = _split_by_sup(G, sup)
            out[f"chr{chrom}"] = ladder_from_panels(site, ref, non)
        except Exception as e:
            out[f"chr{chrom}"] = dict(error=str(e))
    return out


def _split_by_sup(G, sup, per=150, seed=0):
    """Per-superpop 3-way split (matches s16 orientation)."""
    rng = np.random.default_rng(seed)
    si, ri, ni = [], [], []
    for k in np.unique(sup):
        idx = rng.permutation(np.where(sup == k)[0])
        n = idx.size // 3
        a = min(per, n); b = min(per, n); c = min(per, n)
        si.append(idx[:a]); ri.append(idx[a:a + b]); ni.append(idx[a + b:a + b + c])
    take = lambda ii: G[np.concatenate(ii)]
    return take(si), take(ri), take(ni)


def ladder_msprime(seeds=(0, 1, 2)):
    from s18_msprime_pc import simulate_msprime
    acc = {}
    for s in seeds:
        try:
            G, pops = simulate_msprime(seed=s)
        except Exception:
            continue
        site, ref, non = _msprime_split(G, pops)
        r = ladder_from_panels(site, ref, non)
        acc.setdefault("af_only", []).append(r["af_only"])
        acc.setdefault("eig_only", []).append(r["eig_only"])
        for k, v in r["ladder"].items():
            acc.setdefault(k, []).append(v)
    return {str(k): dict(mean=float(np.mean(v)), std=float(np.std(v)))
            for k, v in acc.items()}


def _msprime_split(G, pops, seed=0):
    rng = np.random.default_rng(seed)
    si, ri, ni = [], [], []
    for k in np.unique(pops):
        idx = rng.permutation(np.where(pops == k)[0])
        n = idx.size // 3
        si.append(idx[:n]); ri.append(idx[n:2 * n]); ni.append(idx[2 * n:3 * n])
    take = lambda ii: G[np.concatenate(ii)]
    return take(si), take(ri), take(ni)


PAIRED_CELLS = [(400, 2000, 0.08), (400, 8000, 0.08),
                (1600, 8000, 0.08), (400, 2000, 0.24)]
# the n=1600 cell is BLAS-heavy (SVD inseam ~ n^3); split it out so the fast
# cells aren't blocked behind it (see s29_paired_slow.py).
PAIRED_FAST = [(400, 2000, 0.08), (400, 8000, 0.08)]
PAIRED_SLOW = [(1600, 8000, 0.08), (400, 2000, 0.24)]


def run(shard=0, nshards=1, paired_only=False):
    out_path = os.path.join(OUT, f"s29_negation_shard{shard}.json")
    t0 = time.time()
    res = {}
    # paired sweep: 600 seeds per cell, split across the first 3 shards
    n_sh = min(3, nshards)
    if shard < n_sh:
        per_sh = 600 // n_sh
        seed0 = shard * per_sh
        res["paired"] = paired_sweep(PAIRED_FAST, n_seeds=per_sh, seed0=seed0)
    if not paired_only:
        if shard == 0:
            res["ladder_sim"] = ladder_sim()
        if shard == 1:
            res["ladder_real"] = ladder_real()
        if shard == 2:
            res["ladder_msprime"] = ladder_msprime()
    res["shard"] = shard
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(out_path, "w") as f:
        json.dump(res, f, indent=2)
    print(f"shard{shard} done", res["elapsed_s"], "s", flush=True)


if __name__ == "__main__":
    sh = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    ns = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    po = (len(sys.argv) > 3 and sys.argv[3] == "paired")
    run(sh, ns, paired_only=po)
