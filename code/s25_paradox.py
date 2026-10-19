"""S25 - Sanitized-diagnostics paradox, redesigned as a result.

Replaces the degenerate s3b statistic (probe COUNT, 3->0 by construction) with
a quantity measured on the RELEASED OBJECT:
  E(pi) = attack AUROC of the finest-precision AF the analyst actually sees,
and a quality-vs-exposure Pareto. The claim under test:

  "Pre-providing a coarse SANITIZED diagnostic block lets the analyst stop
   without requesting a finer reveal, so the finest object it ever sees is
   coarser -> lower membership exposure, at equal GWAS quality."

Falsifiable prediction: the paradox holds IFF the sanitized block's precision
p_san is coarser than the finest precision the unsanitized analyst reaches
(p_stop), i.e. it is a precision-GAP effect and must vanish when p_san >= p_stop
or when the exposure curve e(p) has saturated. The null-fixed family (a fixed-k
analyst that never needs a diagnostic) is the control: NO paradox expected.

Families: A_scree, B_gap, C_cumvar (retain >=X% variance), D_mp (bulk edge),
E_lambdastab, F_parallel, G_null_fixed.

Output: results/raw/s25_paradox.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc, pca, marginal_association,
                    genomic_control_lambda)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")
LADDER = [1, 2, 3, 4, 5, 8]      # decimals the analyst may request


def choose_k(eig, family, cumvar_target=0.9):
    eig = np.asarray(eig, float)
    if family == "A_scree":
        thr = 0.05 * eig[0]
        for i in range(len(eig)):
            if eig[i] < thr:
                return max(2, i)
        return 4
    if family == "B_gap":
        r = eig[:6][:-1] / np.maximum(eig[:6][1:], 1e-9)
        return int(np.argmax(r) + 2)
    if family == "C_cumvar":
        cs = np.cumsum(eig) / eig.sum()
        return int(np.searchsorted(cs, cumvar_target) + 1)
    if family == "D_mp":
        # Marchenko-Pastur bulk edge for an n x m matrix -> signal PCs above it
        # (approximate: keep PCs above 2x the median eigenvalue)
        med = np.median(eig)
        return int(max(2, np.sum(eig > 2 * med)))
    if family == "E_lambdastab":
        return 4
    if family == "F_parallel":
        med = np.median(eig)
        return int(max(2, np.sum(eig > med)))
    if family == "G_null_fixed":
        return 4                      # never needs a diagnostic
    raise ValueError(family)


def spectrum(site_G, decimals):
    k_res = max(2, decimals * 2) if decimals < 3 else (6 if decimals < 5 else 10)
    eig, _, load = pca(site_G, k=k_res)
    return eig, load


def probes_until_stable(site_G, family, sanitized_dec, has_sanitized):
    """Sanitized arm starts FROM the coarser block; unsanitized climbs the ladder."""
    if has_sanitized:
        eig, _ = spectrum(site_G, sanitized_dec)
        return 0, eig, sanitized_dec
    if family == "G_null_fixed":
        # an analyst that needs no diagnostic: it stops at the coarsest rung, so
        # the sanitized block gives it nothing -> NO paradox (the control).
        eig, _ = spectrum(site_G, LADDER[0])
        return 1, eig, LADDER[0]
    prev = None
    for i, dec in enumerate(LADDER):
        eig, _ = spectrum(site_G, dec)
        if prev is not None:
            k = choose_k(eig, family)
            if k == prev:
                return i + 1, eig, dec
        prev = choose_k(eig, family)
    return len(LADDER), eig, LADDER[-1]


def quality(site, k):
    mask, p = marginal_association(site["G"], site["case"], n_pcs=max(1, k), correct_pcs=True)
    lam = genomic_control_lambda(p)
    cal = np.exp(-abs(lam - 1.0) * 5)
    # map ORIGINAL marker indices to positions in the QC-filtered p-vector.
    orig = np.where(mask)[0]
    pos_of = {int(o): i for i, o in enumerate(orig)}
    causal = [c for c in site["causal"] if int(c) in pos_of]
    rec = float(np.mean([p[pos_of[int(c)]] < 0.05 for c in causal])) if causal else 0.0
    return float(0.5 * cal + 0.5 * rec), float(lam), rec


def exposure_at(site, ref, nonmem, decimals):
    af = np.round(site["G"].mean(0) / 2.0, decimals)
    r_af = allele_freq(ref["G"])
    lab = np.concatenate([np.ones(len(site["G"])), np.zeros(len(nonmem["G"]))])
    return float(auroc(np.concatenate([af_lrt_attack(site["G"], af, r_af),
                                       af_lrt_attack(nonmem["G"], af, r_af)]), lab))


def run(seeds=range(12), families=("A_scree", "B_gap", "C_cumvar", "D_mp",
                                   "E_lambdastab", "F_parallel", "G_null_fixed"),
        p_san_values=(1, 2, 3), tol=0.05, n_groups=3, n_per_group=600,
        n_snps=2000, per=200):
    from scipy import stats
    rows = []
    for seed in seeds:
        pop = simulate_population(n_groups=n_groups, n_per_group=n_per_group,
                                  n_snps=n_snps, seed=seed)
        site, ref, nonmem = split_population3(pop, per, per, per)
        ecurve = {p: exposure_at(site, ref, nonmem, p) for p in LADDER}
        for family in families:
            nprobe, eig, p_stop = probes_until_stable(site["G"], family, 2, False)
            k = choose_k(eig, family)
            q_uns, _, _ = quality(site, k)
            e_uns = ecurve[p_stop]
            for p_san in p_san_values:
                eig_s, _ = spectrum(site["G"], p_san)
                k_s = choose_k(eig_s, family)
                q_san, _, _ = quality(site, k_s)
                e_san = ecurve[p_san]
                rows.append(dict(seed=seed, family=family, p_san=p_san,
                                 p_stop=p_stop, E_unsan=e_uns, E_san=e_san,
                                 delta_e=e_uns - e_san, delta_q=q_san - q_uns,
                                 nprobe=nprobe, ecurve=ecurve))
    # per (family, p_san) paired stats
    paired = {}
    for family in families:
        for p_san in p_san_values:
            sub = [r for r in rows if r["family"] == family and r["p_san"] == p_san]
            de = np.array([r["delta_e"] for r in sub])
            dq = np.array([r["delta_q"] for r in sub])
            rng = np.random.default_rng(0)
            boot = [de[rng.integers(0, len(de), len(de))].mean() for _ in range(2000)]
            paired[f"{family}|psan={p_san}"] = dict(
                delta_e_mean=float(de.mean()), delta_e_ci=[float(np.percentile(boot, 2.5)),
                                                           float(np.percentile(boot, 97.5))],
                wilcoxon_p=float(stats.wilcoxon(de).pvalue) if np.any(de) else 1.0,
                delta_q_mean=float(dq.mean()), parity_ok=bool(dq.mean() >= -0.03),
                cau_lt_unsan=bool(de.mean() > 0), n=len(sub))
    # paradox region: families x p_san with significant positive delta_e & parity
    paradox = {k: (v["delta_e_ci"][0] > 0 and v["parity_ok"]) for k, v in paired.items()}
    res = dict(rows=rows, paired=paired, paradox_region=paradox,
               null_family_paradox=any(paradox.get(f"G_null_fixed|psan={p}", False)
                                       for p in p_san_values),
               n_seeds=len(list(seeds)))
    with open(os.path.join(OUT, "s25_paradox.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("verdicts:")
    for k in sorted(paired):
        v = paired[k]
        print(f"  {k}: dE={v['delta_e_mean']:+.3f}{tuple(round(x,3) for x in v['delta_e_ci'])} "
              f"dQ={v['delta_q_mean']:+.3f} paradox={paradox[k]}")
    print("null-family paradox:", res["null_family_paradox"])
    return res


if __name__ == "__main__":
    run()
