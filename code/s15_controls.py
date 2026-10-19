"""S15 - supporting controls for numbers cited in the paper's robustness text.

Backs three figures that were previously computed ad hoc:
  (a) the batch-missingness variant vs a full-coverage AF attack and a
      uniform-random-coverage AF attack on the SAME panels (the "0.76 vs 0.845"
      gloss and the coverage-degraded-AF reading);
  (b) the real-vs-sim structure comparison (Fst-like between-group statistic and
      top-10 PC variance fraction) that replaces the retracted "less structured"
      explanation.

Fst-like estimator (named as in the paper): for each SNP, Var_k(p_k) / (p̄(1-p̄))
averaged over SNPs, where p_k is the alt-allele frequency in ancestry group k
and p̄ is the pooled frequency; variance uses numpy's default (ddof=0).

Output: results/raw/s15_controls.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc, pca)
from real_panel import get_real_panel

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def fst_like(G, groups):
    """Between-group variance of allele frequency over the pooled binomial
    variance, averaged over SNPs."""
    af = G.mean(0) / 2.0
    hs = np.mean([2 * (G[groups == g].mean(0) / 2.0) * (1 - G[groups == g].mean(0) / 2.0)
                  for g in np.unique(groups)], axis=0)
    ht = 2 * af * (1 - af)
    num = ht - hs
    return float(np.nanmean(num / np.where(ht > 0, ht, np.nan)))


def top10_frac(G, mask, eig):
    """Fraction of TOTAL genetic variance (all PCs) captured by the top-10 PCs."""
    X = G[:, mask].astype(float); X = X - X.mean(0)
    S = np.linalg.svd(X, compute_uv=False)
    total_eig = (S ** 2).sum() / max(1, X.shape[0] - 1)
    return float(eig.sum() / total_eig)


def batch_controls(seeds=(0, 2, 4)):  # reuse the S9a seed set
    rows = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, non = split_population3(pop, 300, 300, 300)
        mask = release_objects(site)["mask"]
        af_s, af_r = allele_freq(site["G"]), allele_freq(ref["G"])
        y = np.concatenate([np.ones(len(site["G"])), np.zeros(len(non["G"]))])

        def auc(sm, sn):
            return float(auroc(np.concatenate([sm, sn]), y))

        full = auc(af_lrt_attack(site["G"], af_s, af_r, mask=mask),
                   af_lrt_attack(non["G"], af_s, af_r, mask=mask))
        # uniform random 70% coverage of SNPs (shared by members and non-members)
        rng = np.random.default_rng(seed + 101)
        cov = rng.random(site["G"].shape[1]) < 0.70
        vis = mask & cov
        unif = auc(af_lrt_attack(site["G"], af_s, af_r, mask=vis),
                   af_lrt_attack(non["G"], af_s, af_r, mask=vis))
        # batch-mixed: members/non-members split across two marker batches
        batchA = rng.random(site["G"].shape[1]) < 0.70
        mem_b = rng.random(site["G"].shape[0]) < 0.5
        non_b = rng.random(non["G"].shape[0]) < 0.5

        def score(panel, flag):
            out = np.zeros(panel["G"].shape[0])
            for i in range(panel["G"].shape[0]):
                v = mask & (batchA if flag[i] else ~batchA)
                out[i] = float(af_lrt_attack(panel["G"][i], af_s, af_r, mask=v)[0])
            return out
        mixed = float(auroc(np.concatenate([score(site, mem_b), score(non, non_b)]), y))
        rows.append(dict(seed=seed, full_coverage_af=full,
                         uniform70_coverage_af=unif, batch_mixed=mixed))
        print("batch", json.dumps(rows[-1]), flush=True)
    return rows


def structure_compare():
    # real
    G, sup, pos = get_real_panel()
    rng = np.random.default_rng(0)
    idx = rng.permutation(np.arange(G.shape[0]))
    per = G.shape[0] // 3
    site = G[idx[:per]]; ss = sup[idx[:per]]
    mask = np.ones(G.shape[1], bool)
    eig_r = pca(site, k=10)[0]
    real = dict(fst_like=fst_like(site, ss), top10_var_frac=top10_frac(site, mask, eig_r))
    # sim matched (450/group, 1363 SNPs, ld_blocks=100)
    pop = simulate_population(n_groups=3, n_per_group=450, n_snps=1363, ld_blocks=100, seed=0)
    s, r, n = split_population3(pop, 150, 150, 150)
    rel = release_objects(s)
    sim = dict(fst_like=fst_like(s["G"], s["ancestry"]),
               top10_var_frac=top10_frac(s["G"], rel["mask"], rel["eig"]))
    print("structure real", real, "sim", sim, flush=True)
    return dict(real=real, sim_matched=sim)


def run():
    res = dict(batch_controls=dict(per_seed=batch_controls()), structure=structure_compare())
    bc = res["batch_controls"]["per_seed"]
    res["batch_controls"]["agg"] = {
        k: float(np.mean([r[k] for r in bc])) for k in
        ("full_coverage_af", "uniform70_coverage_af", "batch_mixed")}
    with open(os.path.join(OUT, "s15_controls.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("agg", json.dumps({k: round(v, 4) for k, v in res["batch_controls"]["agg"].items()}),
          flush=True)
    return res


if __name__ == "__main__":
    run()
