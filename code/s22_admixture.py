"""S22 - admixture stress test.

The paper's simulator (and the msprime cross-check) uses DISCRETE population
structure. Real biobank cohorts are admixed. This tests whether PC dominance
requires discrete structure: under continuous admixture the admixture gradient
is a top PC shared equally by members and non-members, so the in-sample residual
share should shrink and the PC channel should attenuate toward the real-data
regime while the AF-LRT channel (a finite-sample study-vs-ref effect) stays
roughly flat.

Two arms:
  - parametric: per-individual Dirichlet ancestry weights w_i ~ Dirichlet(alpha)
    over 3 source groups, f_i = sum_k w_ik * group_freq_k, G_i ~ Binomial(2, f_i).
  - msprime: pulse admixture (ADC = 0.5*A + 0.5*B) with real recombination.

Output: results/raw/s22_admixture.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import allele_freq, af_lrt_attack, auroc, pca

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def simulate_admixed(n=900, n_snps=4000, fst=0.08, alpha=1.0, seed=0, n_src=3):
    rng = np.random.default_rng(seed)
    base = rng.uniform(0.05, 0.95, size=n_snps)
    if fst <= 0:
        gf = np.tile(base, (n_src, 1))
    else:
        a = base * (1 - fst) / fst; b = (1 - base) * (1 - fst) / fst
        gf = rng.beta(a[None, :], b[None, :], size=(n_src, n_snps))
    W = rng.dirichlet(np.full(n_src, alpha), size=n)         # (n, K)
    F = W @ gf                                                # (n, m) per-individual freq
    G = rng.binomial(2, np.clip(F, 1e-4, 1 - 1e-4)).astype(np.int8)
    # LD injection (adjacent) so the panels resemble the parametric regime
    for b0 in range(0, n_snps - 6, 6):
        latent = rng.normal(0, 1, size=n)
        hi = latent > 0.5
        for c in range(b0, min(b0 + 6, n_snps)):
            flip = rng.random(n) < 0.15
            G[hi & flip, c] = np.minimum(2, G[hi & flip, c] + 1)
    dom = W.argmax(1)     # modal ancestry, for the disjoint non-member split
    return dict(G=G, ancestry=dom, W=W)


def split_three(pop, per):
    G, anc = pop["G"], pop["ancestry"]
    sidx, ridx, nidx = [], [], []
    for k in np.unique(anc):
        idx = np.where(anc == k)[0]
        q = min(per, idx.size // 3)
        sidx.append(idx[:q]); ridx.append(idx[q:2 * q]); nidx.append(idx[2 * q:3 * q])
    sub = lambda i: G[np.concatenate(i)]
    return sub(sidx), sub(ridx), sub(nidx)


def pc_dist(G, mask, L, mean, k=None):
    Lk = L if k is None else L[:k]
    X = G[:, mask].astype(float) - mean[None, :]
    return np.linalg.norm((X @ Lk[:, mask].T) / np.sqrt(max(1, mask.sum())), axis=1)


def cal(a):
    return max(float(a), 1 - float(a))


def evaluate(site, ref, non, seed):
    mask = np.ones(site.shape[1], bool)
    _, _, L = pca(site, k=10); mean = site.mean(0)
    y = np.concatenate([np.ones(len(site)), np.zeros(len(non))])
    af_s, af_r = allele_freq(site), allele_freq(ref)
    af = float(auroc(np.concatenate([af_lrt_attack(site, af_s, af_r, mask=mask),
                                     af_lrt_attack(non, af_s, af_r, mask=mask)]), y))
    dm, dn = pc_dist(site, mask, L, mean), pc_dist(non, mask, L, mean)
    pc = cal(auroc(np.concatenate([dm, dn]), y))
    for K in (1, 3, 5):
        pass
    top1 = cal(auroc(np.concatenate([pc_dist(site, mask, L, mean, 1),
                                     pc_dist(non, mask, L, mean, 1)]), y))
    top5 = cal(auroc(np.concatenate([pc_dist(site, mask, L, mean, 5),
                                     pc_dist(non, mask, L, mean, 5)]), y))
    rng = np.random.default_rng(seed)
    Lr = rng.normal(size=L.shape); Lr, _ = np.linalg.qr(Lr.T); Lr = Lr.T
    eig = cal(auroc(np.concatenate([pc_dist(site, mask, Lr, mean),
                                    pc_dist(non, mask, Lr, mean)]), y))
    return dict(af_only=af, pc_calibrated=pc, top1=top1, top5=top5, eig_only=eig,
                delta_pc_minus_af=pc - af)


def parametric_arm(alphas=(0.2, 1.0, 10.0), n_seeds=10, per=200):
    out = {}
    for alpha in alphas:
        rows = []
        for seed in range(n_seeds):
            pop = simulate_admixed(n=3 * per * 3, alpha=alpha, seed=seed)
            site, ref, non = split_three(pop, per)
            rows.append(evaluate(site, ref, non, seed))
        agg = {k: [float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows]))]
               for k in rows[0]}
        out[f"alpha_{alpha}"] = dict(agg=agg, per_seed=rows)
        print("param", alpha, json.dumps({k: round(v[0], 4) for k, v in agg.items()}), flush=True)
    return out


def msprime_arm(seeds=(0, 1, 2), n_per=200, seq_len=5_000_000, n_snps=8000):
    import msprime
    rows = []
    for seed in seeds:
        dem = msprime.Demography()
        for name in ("A", "B", "ADC"):
            dem.add_population(name=name, initial_size=10_000)
        dem.add_population(name="ANC", initial_size=10_000)
        dem.add_population_split(time=2000, derived=["A", "B"], ancestral="ANC")
        dem.set_symmetric_migration_rate(["A", "B"], 1e-4)
        dem.add_admixture(time=100, derived="ADC", ancestral=["A", "B"],
                          proportions=[0.5, 0.5])
        dem.sort_events()
        ts = msprime.sim_ancestry(samples={"A": n_per, "B": n_per, "ADC": n_per},
                                  demography=dem, sequence_length=seq_len,
                                  recombination_rate=1.25e-8, random_seed=seed + 1)
        ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=seed + 2)
        G = ts.genotype_matrix().T.astype(np.int8)
        pops = np.array([ts.tables.nodes[i].population if ts.tables.nodes[i].population >= 0
                         else 0 for i in ts.samples()])
        maf = np.minimum(G.mean(0) / 2.0, 1 - G.mean(0) / 2.0)
        keep = np.where(maf > 0.05)[0]
        rng = np.random.default_rng(seed)
        if keep.size > n_snps:
            keep = np.sort(rng.choice(keep, size=n_snps, replace=False))
        G = G[:, keep]; pops = pops
        # DNA of the "study" = admixed ADC + one source; non-members = the other
        # Use within-population disjoint split of the whole admixed panel.
        rng2 = np.random.default_rng(seed + 9)
        idx = rng2.permutation(np.arange(G.shape[0]))
        q = G.shape[0] // 3
        site, ref, non = G[idx[:q]], G[idx[q:2 * q]], G[idx[2 * q:3 * q]]
        rows.append(evaluate(site, ref, non, seed))
        print("msprime-admix", seed, json.dumps({k: round(v, 4) for k, v in rows[-1].items()}),
              flush=True)
    agg = {k: [float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows]))]
           for k in rows[0]}
    return dict(agg=agg, per_seed=rows)


def run():
    res = dict(parametric=parametric_arm(), msprime_admix=msprime_arm())
    with open(os.path.join(OUT, "s22_admixture.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
