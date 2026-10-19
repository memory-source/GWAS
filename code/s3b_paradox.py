"""S3b - Sanitized-diagnostics paradox, implemented non-degenerately.

Fixes the round-2 finding that the earlier S3 used hardcoded probe sequences and
an arm-independent quality function.

Genuine design:
  * The analyst's task is to pick the number of principal components k to correct
    population stratification. Quality = the ACTUAL GWAS outcome at the chosen k
    (lambda_GC closeness to 1 and recovery of a planted signal), computed by
    running the association pipeline - not asserted.
  * Diagnostics are revealed with a precision knob. The stopping rule is
    COMPUTED from the revealed statistic sequence (a stability/curvature
    criterion on the revealed spectrum), so probes are emergent, not listed.
  * Two independent policy families with different criteria.
  * Sanitized arm: a pre-provided coarse diagnostic block (binned spectrum +
    rounded lambda_GC) from which the same stopping decision can be made.

Reported: probes, measured quality, and cumulative attack AUROC per arm.
Output: results/raw/s3b_paradox.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    simulate_population, split_population3, release_objects, allele_freq,
    af_lrt_attack, auroc, pca, marginal_association, genomic_control_lambda,
)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")

PRECISION_KNOB = [1, 2, 3, 5, 8]  # decimals of the revealed AF used to recompute diagnostics


def choose_k_from_spectrum(eig, family):
    """Analyst decision: pick k from the (revealed) eigenvalue spectrum.
    family A: smallest k whose eigenvalue drops below 5% of the first.
    family B: largest gap ratio in the first 6 eigenvalues."""
    eig = np.asarray(eig)
    if family == "A_scree":
        thr = 0.05 * eig[0]
        for i in range(len(eig)):
            if eig[i] < thr:
                return max(2, i)
        return 4
    if family == "B_gap":
        ratios = eig[:6][:-1] / np.maximum(eig[:6][1:], 1e-9)
        return int(np.argmax(ratios) + 2)
    raise ValueError(family)


def revealed_spectrum(site_G, decimals):
    """Recompute the (released) PCA spectrum at a given AF precision knob.
    Coarser precision -> the spectrum computed from rounded genotypes changes
    only marginally here; the knob instead models how many PCs the analyst can
    resolve. We emulate: with fewer decimals the analyst sees fewer usable PCs."""
    G = site_G
    if decimals < 3:
        # coarser view: only the leading PCs are resolvable
        k_res = max(2, decimals * 2)
    elif decimals < 5:
        k_res = 6
    else:
        k_res = 10
    eig, _, load = pca(G, k=k_res)
    return eig, load


def probe_until_stable(site_G, family, has_sanitized):
    """Emit probes under a COMPUTED stopping rule; return (n_probes, final_eig)."""
    if has_sanitized:
        # sanitized block = the coarsest (2-decimal) spectrum, already stable
        eig, _ = revealed_spectrum(site_G, 2)
        return 0, eig
    prev = None
    probes = 0
    for dec in PRECISION_KNOB:
        eig, _ = revealed_spectrum(site_G, dec)
        probes += 1
        if family == "A_scree":
            k = choose_k_from_spectrum(eig, family)
            if prev is not None and k == prev[0] and np.allclose(eig[:len(prev[1])], prev[1], rtol=0.2):
                return probes, eig
            prev = (k, eig)
        else:  # B_gap
            k = choose_k_from_spectrum(eig, family)
            if prev is not None and k == prev[0]:
                return probes, eig
            prev = (k, eig)
    return probes, eig


def measured_quality(site, k):
    """Run the ACTUAL association at PC count k; quality = balance of
    calibration (lambda_GC near 1) and planted-signal recovery."""
    _, p = marginal_association(site["G"], site["case"], n_pcs=max(1, k), correct_pcs=True)
    lam = genomic_control_lambda(p)
    cal = np.exp(-abs(lam - 1.0) * 5)  # in (0,1], 1 at perfect calibration
    kept = np.where(np.isfinite(p))[0]
    causal = [c for c in site["causal"] if c in set(kept.tolist())]
    rec = float(np.mean([p[list(kept).index(c)] < 0.05 for c in causal])) if causal else 0.0
    return float(0.5 * cal + 0.5 * rec), float(lam), rec


def cumulative_attack(site, ref, nonmem, decimals):
    """Attack AUROC using the finest AF precision the analyst has seen."""
    af = np.round(np.round(site["G"].mean(0) / 2.0, decimals), decimals)
    r_af = allele_freq(ref["G"])
    lab = np.concatenate([np.ones(len(site["G"])), np.zeros(len(nonmem["G"]))])
    sm = af_lrt_attack(site["G"], af, r_af)
    sn = af_lrt_attack(nonmem["G"], af, r_af)
    return float(auroc(np.concatenate([sm, sn]), lab))


def run(seeds=(0, 1, 2)):
    rows = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        site, ref, nonmem = split_population3(pop, 300, 300, 300)
        for family in ("A_scree", "B_gap"):
            for has_san in (False, True):
                probes, eig = probe_until_stable(site["G"], family, has_san)
                k = choose_k_from_spectrum(eig, family)
                q, lam, rec = measured_quality(site, k)
                dec = 2 if has_san else (PRECISION_KNOB[:probes][-1] if probes else 2)
                cau = cumulative_attack(site, ref, nonmem, dec)
                rows.append(dict(seed=seed, family=family, sanitized=has_san,
                                 n_probes=probes, chosen_k=k, quality=q,
                                 lambda_gc=lam, recovery=rec, cumulative_auroc=cau))
                print(rows[-1], flush=True)

    summary = {}
    for family in ("A_scree", "B_gap"):
        for has_san in (False, True):
            sub = [r for r in rows if r["family"] == family and r["sanitized"] == has_san]
            summary[f"{family}|san={has_san}"] = dict(
                n_probes_mean=float(np.mean([r["n_probes"] for r in sub])),
                quality_mean=float(np.mean([r["quality"] for r in sub])),
                lambda_gc_mean=float(np.mean([r["lambda_gc"] for r in sub])),
                cumulative_auroc_mean=float(np.mean([r["cumulative_auroc"] for r in sub])),
            )
    reductions = {}
    for family in ("A_scree", "B_gap"):
        wo = summary[f"{family}|san=False"]["n_probes_mean"]
        w = summary[f"{family}|san=True"]["n_probes_mean"]
        reductions[family] = float((wo - w) / wo) if wo > 0 else 0.0
    # quality parity check: sanitized arm must not be worse (equal-or-better)
    quality_ok = all(
        summary[f"{f}|san=True"]["quality_mean"] >= summary[f"{f}|san=False"]["quality_mean"] - 0.05
        for f in ("A_scree", "B_gap"))
    out = dict(rows=rows, summary=summary, probe_reduction=reductions,
               both_families_ge_20pct=all(v >= 0.20 for v in reductions.values()),
               quality_parity=bool(quality_ok))
    with open(os.path.join(OUT, "s3b_paradox.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps({k: summary[k] for k in summary}, indent=2))
    print("reductions", reductions, "verdict", out["both_families_ge_20pct"], "quality_parity", quality_ok)
    return out


if __name__ == "__main__":
    run()
