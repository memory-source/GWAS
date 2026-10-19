"""S4 - Realistic-population cross-check with msprime (physics-based simulator).

Builds a 3-population structured demography with migration over a 2 Mb region
(real recombination, realistic site-frequency spectrum and LD), extracts
genotypes, and reruns the core leakage measurements. Confirms the parametric
simulator's trend direction on an independent, more realistic data source.

Output: results/raw/s4_msprime.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    release_objects, allele_freq, af_lrt_attack, auroc, tpr_at_fpr,
)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")
os.makedirs(OUT, exist_ok=True)


def simulate_msprime(n_per_pop=300, seq_len=2_000_000, seed=0, n_snps=4000):
    import msprime
    dem = msprime.Demography()
    dem.add_population(name="A", initial_size=10_000)
    dem.add_population(name="B", initial_size=10_000)
    dem.add_population(name="C", initial_size=10_000)
    dem.add_population(name="ANC", initial_size=10_000)
    dem.add_population_split(time=2000, derived=["A", "B", "C"], ancestral="ANC")
    dem.set_symmetric_migration_rate(["A", "B", "C"], 1e-4)
    ts = msprime.sim_ancestry(
        samples={"A": n_per_pop, "B": n_per_pop, "C": n_per_pop},
        demography=dem, sequence_length=seq_len, recombination_rate=1.25e-8,
        random_seed=seed + 1,
    )
    ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=seed + 2)
    G = ts.genotype_matrix().T  # (n_individuals, n_variants)
    pops = np.array([ts.tables.nodes[i].population if ts.tables.nodes[i].population >= 0
                     else 0 for i in ts.samples()])
    G = G.astype(np.int8)
    # downsample variants to n_snps by MAF filtering + random selection
    # (G.mean(0) is the ALT allele frequency times 2)
    maf = np.minimum(G.mean(0) / 2.0, 1 - G.mean(0) / 2.0)
    keep = np.where((maf > 0.05))[0]
    rng = np.random.default_rng(seed)
    if keep.size > n_snps:
        keep = rng.choice(keep, size=n_snps, replace=False)
        keep.sort()
    G = G[:, keep]
    return G, pops, ts.num_sites


def run(seeds=(0, 1, 2)):
    res = {"per_seed": []}
    for seed in seeds:
        G, pops, n_sites = simulate_msprime(seed=seed)
        # three-way disjoint split: study / attacker reference / non-member targets
        sidx, ridx, nidx = [], [], []
        for k in np.unique(pops):
            idx = np.where(pops == k)[0]
            sidx.append(idx[:100])
            ridx.append(idx[100:200])
            nidx.append(idx[200:300])
        sidx = np.concatenate(sidx)
        ridx = np.concatenate(ridx)
        nidx = np.concatenate(nidx)
        study = dict(G=G[sidx], ancestry=pops[sidx], case=np.zeros(len(sidx), dtype=np.int8),
                     causal=np.array([], dtype=int))
        ref = dict(G=G[ridx], ancestry=pops[ridx], case=np.zeros(len(ridx), dtype=np.int8),
                   causal=np.array([], dtype=int))
        nonmem = dict(G=G[nidx], ancestry=pops[nidx], case=np.zeros(len(nidx), dtype=np.int8),
                      causal=np.array([], dtype=int))
        maf = np.minimum(G.mean(0) / 2.0, 1 - G.mean(0) / 2.0)
        mask = (maf > 0.05) & (maf < 0.95)
        af_s = allele_freq(study["G"])
        af_r = allele_freq(ref["G"])

        sm = af_lrt_attack(study["G"], af_s, af_r, mask=mask)
        sn = af_lrt_attack(nonmem["G"], af_s, af_r, mask=mask)
        labels = np.concatenate([np.ones(len(sm)), np.zeros(len(sn))])
        sc = np.concatenate([sm, sn])
        row = dict(seed=seed, n_snps=int(mask.sum()), n_sites=int(n_sites),
                   auroc_full=float(auroc(sc, labels)),
                   tpr_at_fpr1=float(tpr_at_fpr(sc, labels, 0.01)),
                   tpr_at_fpr01=float(tpr_at_fpr(sc, labels, 0.001)))
        # rounded (F4-style) precision, on the INDEPENDENT non-member panel
        smr = af_lrt_attack(study["G"], np.round(af_s, 1), af_r, mask=mask)
        snr = af_lrt_attack(nonmem["G"], np.round(af_s, 1), af_r, mask=mask)
        scr = np.concatenate([smr, snr])
        row["auroc_round1"] = float(auroc(scr, labels))
        row["tpr_at_fpr1_round1"] = float(tpr_at_fpr(scr, labels, 0.01))
        res["per_seed"].append(row)
        print(json.dumps(row))

    for k in ("auroc_full", "auroc_round1", "tpr_at_fpr1", "tpr_at_fpr1_round1"):
        res[k + "_mean"] = float(np.mean([r[k] for r in res["per_seed"]]))
        res[k + "_std"] = float(np.std([r[k] for r in res["per_seed"]]))
    with open(os.path.join(OUT, "s4_msprime.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
