"""S18 - msprime cross-check for the PC channel (realistic recombination/LD).

The PC-channel experiments all used the parametric simulator with injected LD
blocks. This closes that gap: build an msprime 3-population demography with real
recombination over a long region, extract genotypes, and run the same
orientation + mechanism + granularity measurements, so the PC channel is
characterized on a physics-based simulator as well.

Output: results/raw/s18_msprime_pc.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import allele_freq, af_lrt_attack, auroc, pca

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def cal(a):
    return max(float(a), 1 - float(a))


def pc_dist(G, loadings, mean, k=None):
    L = loadings if k is None else loadings[:k]
    X = G.astype(np.float64) - mean[None, :]
    return np.linalg.norm((X @ L.T) / np.sqrt(max(1, G.shape[1])), axis=1)


def simulate_msprime(n_per_pop=400, seq_len=5_000_000, seed=0, n_snps=8000):
    import msprime
    dem = msprime.Demography()
    for name in ("A", "B", "C"):
        dem.add_population(name=name, initial_size=10_000)
    dem.add_population(name="ANC", initial_size=10_000)
    dem.add_population_split(time=2000, derived=["A", "B", "C"], ancestral="ANC")
    dem.set_symmetric_migration_rate(["A", "B", "C"], 1e-4)
    ts = msprime.sim_ancestry(
        samples={"A": n_per_pop, "B": n_per_pop, "C": n_per_pop},
        demography=dem, sequence_length=seq_len, recombination_rate=1.25e-8,
        random_seed=seed + 1)
    ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=seed + 2)
    G = ts.genotype_matrix().T.astype(np.int8)
    pops = np.array([ts.tables.nodes[i].population if ts.tables.nodes[i].population >= 0
                     else 0 for i in ts.samples()])
    maf = np.minimum(G.mean(0) / 2.0, 1 - G.mean(0) / 2.0)
    keep = np.where(maf > 0.05)[0]
    rng = np.random.default_rng(seed)
    if keep.size > n_snps:
        keep = np.sort(rng.choice(keep, size=n_snps, replace=False))
    return G[:, keep], pops


def three_way(G, pops, seed):
    sidx, ridx, nidx = [], [], []
    for k in np.unique(pops):
        idx = np.where(pops == k)[0]
        n = idx.size // 3
        sidx.append(idx[:n]); ridx.append(idx[n:2 * n]); nidx.append(idx[2 * n:3 * n])
    return G[np.concatenate(sidx)], G[np.concatenate(ridx)], G[np.concatenate(nidx)]


def studies(G, pops, seed):
    site, ref, non = three_way(G, pops, seed)
    mask = np.ones(G.shape[1], bool)
    _, _, L = pca(site, k=10); mean = site.mean(0)
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    af_s, af_r = allele_freq(site), allele_freq(ref)
    am = af_lrt_attack(site, af_s, af_r, mask=mask)
    an = af_lrt_attack(non, af_s, af_r, mask=mask)
    dm, dn = pc_dist(site, L, mean), pc_dist(non, L, mean)
    out = dict(seed=seed, af_only=float(auroc(np.concatenate([am, an]), y)),
               pc_naive=float(auroc(np.concatenate([dm, dn]), y)),
               pc_calibrated=cal(auroc(np.concatenate([dm, dn]), y)))
    # granularity
    for K in (1, 3, 5):
        dmK, dnK = pc_dist(site, L, mean, k=K), pc_dist(non, L, mean, k=K)
        out[f"top{K}"] = cal(auroc(np.concatenate([dmK, dnK]), y))
    out["full"] = out["pc_calibrated"]
    rng = np.random.default_rng(seed)
    Lr = rng.normal(size=L.shape); Lr, _ = np.linalg.qr(Lr.T); Lr = Lr.T
    out["eig_only"] = cal(auroc(np.concatenate([pc_dist(site, Lr, mean),
                                                pc_dist(non, Lr, mean)]), y))
    # held-out-site mechanism: WITHIN-population disjoint split so members,
    # the other site, and the shared non-members are ancestry-matched (the S9d /
    # S14 design). Splitting by msprime population instead would confound the
    # null with ancestry separation.
    rng = np.random.default_rng(seed)
    A_i, B_i, C_i = [], [], []
    for k in np.unique(pops):
        idx = rng.permutation(np.where(pops == k)[0])
        q = idx.size // 4
        A_i.append(idx[:q]); B_i.append(idx[q:2 * q]); C_i.append(idx[2 * q:3 * q])
    A, B, C = G[np.concatenate(A_i)], G[np.concatenate(B_i)], G[np.concatenate(C_i)]
    _, _, LA = pca(A, k=10); mA = A.mean(0)
    _, _, LB = pca(B, k=10); mB = B.mean(0)
    yc = np.concatenate([np.ones(len(A)), np.zeros(len(C))])
    out["A_in"] = cal(auroc(np.concatenate([pc_dist(A, LA, mA), pc_dist(C, LA, mA)]), yc))
    out["A_oos"] = cal(auroc(np.concatenate([pc_dist(A, LB, mB), pc_dist(C, LB, mB)]), yc))
    out["A_swapped"] = cal(auroc(np.concatenate([pc_dist(B, LA, mA), pc_dist(C, LA, mA)]), yc))
    return out


def run(seeds=(0, 1, 2, 3, 4)):
    rows = []
    for s in seeds:
        G, pops = simulate_msprime(seed=s)
        r = studies(G, pops, s)
        rows.append(r)
        print(json.dumps({k: round(v, 4) if isinstance(v, float) else v for k, v in r.items()}),
              flush=True)
    keys = ["af_only", "pc_calibrated", "top1", "top3", "top5", "full", "eig_only",
            "A_in", "A_oos", "A_swapped"]
    agg = {k: (float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows])))
           for k in keys}
    res = dict(per_seed=rows, agg=agg, n_seeds=len(seeds))
    with open(os.path.join(OUT, "s18_msprime_pc.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("agg", json.dumps({k: round(v[0], 4) for k, v in agg.items()}), flush=True)
    return res


if __name__ == "__main__":
    run()
