"""S71 - scientific cost under protection (round 22).

The GPT analysis asks that the honest GWAS be characterized beyond lambda-GC and the
minimum p-value: false-positive rate, power, effect agreement, individual QC, under
different release strategies.

Release strategies compared on the honest study's own association analysis:
  - baseline  : the honest study analyzed alone (its own top-10 PCs), no partner release;
  - passive   : the partner input never changes; honest analysis re-fit each round with
                the joint-PCA covariates of (honest, partner);
  - perturbed : the partner's AF-preserving block swap each round; honest analysis
                re-fit with the per-round joint-PCA covariates.

The phenotype is simulated with EXPOSED causal markers and EXPOSED per-marker effects
(eff ~ N(0, 0.35)); the effect comparison is the scale-free Pearson correlation between
the estimated per-marker effect and the true effect ACROSS THE CAUSAL MARKERS, so it is
an effect-agreement measure, not a bias in a fixed unit.

For each arm the script reports, over 8 rounds and 5 seeds:
  - power    = fraction of the planted causal markers recovered at p < 1e-4;
  - off_causal_fraction = fraction of the tested off-causal markers below the same
    threshold; this is a discovery count on markers that were NOT planted, NOT a clean
    null false-positive rate: under linkage disequilibrium some off-causal markers are
    proxies of causal ones and their discoveries are not false;
  - lambda_GC;
  - effect_agreement = Pearson r(estimated effect, true effect) over causal markers;
  - individual QC = max abs change in the partner's per-SNP variance and per-row
    heterozygosity over the release sequence (histograms identical by construction).

Output: results/raw/s71_protection_cost_chr<C>.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
import s51_real_lowstruct as s51
import s61_spectral_tomography as s61
from common import genomic_control_lambda, qc_mask, pca
from real_multi import SUPERPOP


def _pca_scores(G, k):
    """Per-individual top-k PC scores of G (the baseline covariates)."""
    _, scores, _ = pca(G, k=k)
    return scores

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
N_ROUNDS = 8
P_THRESH = 1e-4
N_CAUSAL = 8


def cohorts_for(chrom, ld_window=1500, max_snps=20000):
    Ge, _, _ = s48.load_for(("EUR",), chrom, ld_window, max_snps)
    G3, _, sup3, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    return {"within:EUR": (Ge, np.array(["EUR"] * len(Ge))),
            "super:3pop": (G3, sup3)}


def sim_phenotype_exposed(site, seed, n_causal=N_CAUSAL):
    """Case/control with the causal markers and their per-marker effects returned."""
    rng = np.random.default_rng(seed + 900)
    m = site.shape[1]
    causal = rng.choice(m, size=n_causal, replace=False)
    eff = rng.normal(0.0, 0.35, size=n_causal)
    risk = site[:, causal].astype(float) @ eff
    risk = risk - risk.mean()
    prob = 1.0 / (1.0 + np.exp(-risk))
    case = (rng.random(site.shape[0]) < prob).astype(np.int8)
    return case, causal, eff


def _scan(honest, case, scores, k, mask):
    """GWAS scan on the honest rows with the supplied covariate scores (n x k),
    restricted to the QC-passing marker set `mask`. Returns FULL-LENGTH (M,) arrays
    `pv` and `beta`, with the non-tested markers set to 1.0 / 0.0, so every marker
    index aligns with `causal`."""
    from scipy.stats import chi2 as _chi2
    m = honest.shape[1]
    design = np.column_stack([np.ones(len(honest))] + [scores[:, j] for j in range(k)])
    X = honest[:, mask].astype(float)
    bx, *_ = np.linalg.lstsq(design, X, rcond=None)
    Xr = X - design @ bx
    y = case.astype(float)
    by, *_ = np.linalg.lstsq(design, y, rcond=None)
    yr = y - design @ by
    xc = Xr - Xr.mean(0); yc = yr - yr.mean()
    xss = (xc ** 2).sum(0)
    den = np.sqrt(xss * (yc ** 2).sum()); den[den == 0] = np.inf
    r = np.clip((xc * yc[:, None]).sum(0) / den, -1 + 1e-12, 1 - 1e-12)
    pv_kept = _chi2.sf(r ** 2 * len(honest), 1)
    beta_kept = np.where(xss > 0, (xc * yc[:, None]).sum(0) / np.where(xss > 0, xss, 1.0), 0.0)
    pv = np.ones(m); beta = np.zeros(m)
    pv[mask] = pv_kept; beta[mask] = beta_kept
    return pv, beta


def _summarize(pv, beta, causal, eff, mask):
    """Summarize a full-length scan against the causal markers, using ONLY the tested
    markers (`mask`) for the calibration statistics (lambda-GC over the tested
    markers; off-causal discovery denominator = tested off-causal markers; effect
    agreement over tested causal markers). Power counts ALL causal markers (a
    QC-excluded causal marker is a miss), so it is comparable to a real GWAS."""
    m = len(pv)
    tested = np.asarray(mask, bool)
    hit = np.zeros(m, bool); hit[causal] = True
    recovered = (pv < P_THRESH) & hit & tested
    fp = (pv < P_THRESH) & tested & (~hit)
    n_causal = int(hit.sum())
    n_test_non = int((tested & ~hit).sum())
    causal_tested = tested[causal]                       # length len(causal)
    agreement = (float(np.corrcoef(beta[causal][causal_tested], eff[causal_tested])[0, 1])
                 if causal_tested.sum() > 2 else float("nan"))
    minp = (float(np.min(pv[causal][causal_tested])) if causal_tested.any() else float("nan"))
    return {"lambda_gc": float(genomic_control_lambda(pv[tested])),
            "power": float(recovered.sum()) / n_causal if n_causal else float("nan"),
            "n_causal_recovered": int(recovered.sum()), "n_causal": n_causal,
            "n_causal_in_qc_mask": int(causal_tested.sum()),
            "off_causal_fraction": float(fp.sum()) / max(1, n_test_non),
            "n_off_causal_discoveries": int(fp.sum()), "n_test_off_causal": n_test_non,
            "effect_agreement": agreement,
            "min_p_causal": minp}


def one_split(G, groups, seed, k=s61.K_TOP, block_len=s61.BLOCK_LEN, caps=None):
    roles = s61.split_roles(groups, caps or s61.CAPS, seed)
    honest, partner = G[roles["honest"]], G[roles["partner"]]
    mu, scale = s61.joint_baseline(honest, partner)
    Zh = s61.zscore(honest, mu, scale)
    n_blocks = G.shape[1] // block_len
    case, causal, eff = sim_phenotype_exposed(honest, seed)
    # ONE marker universe for every arm (the QC mask), so all arms are comparable
    mask = qc_mask(honest)

    # baseline: honest alone with its own top-10 PCs (same mask)
    _, _, Ut = None, None, None
    scores0 = _pca_scores(honest, k)
    pv_b, beta_b = _scan(honest, case, scores0, k, mask)
    base = _summarize(pv_b, beta_b, causal, eff, mask)

    out = {"seed": seed, "baseline": base, "n_markers_tested": int(mask.sum()),
           "rounds": {"passive": [], "perturbed": []}, "qc": {}}
    sig0 = np.sort(partner, axis=0)
    var0 = partner.astype(float).var(0)
    het0 = (partner == 1).mean(1)
    max_het = 0.0; max_var = 0.0; hist_ok = True
    rng = np.random.default_rng(seed + 101)
    for strat in ("passive", "perturbed"):
        for t in range(N_ROUNDS):
            if strat == "passive" or t == 0:
                own = partner
            else:
                own, _ = s61.swap_block(partner, rng.choice(partner.shape[0], 2,
                                                            replace=False),
                                        int(rng.integers(0, n_blocks)), block_len)
            hist_ok = hist_ok and bool(np.array_equal(np.sort(own, axis=0), sig0))
            max_var = max(max_var, float(np.max(np.abs(own.astype(float).var(0) - var0))))
            max_het = max(max_het, float(np.max(np.abs((own == 1).mean(1) - het0))))
            Zp = s61.zscore(own, mu, scale)
            lam, V = s61.topk([Zh, Zp], k)
            scores = Zh @ V
            pv, beta = _scan(honest, case, scores, k, mask)
            out["rounds"][strat].append(_summarize(pv, beta, causal, eff, mask))
    out["qc"] = {"histograms_identical": bool(hist_ok),
                 "max_abs_variance_change": max_var,
                 "max_abs_row_heterozygosity_change": max_het}
    return out


def _mean(path, key):
    vals = [r[key] for r in path]
    vals = [v for v in vals if v is not None and np.isfinite(v)]
    return float(np.mean(vals)) if vals else float("nan")


def run(chrom=22, seeds=5):
    out = {"chrom": chrom, "seeds": seeds, "p_threshold": P_THRESH, "cohorts": {}}
    for name, (G, groups) in cohorts_for(chrom).items():
        rows = [one_split(G, groups, s) for s in range(seeds)]
        blk = {"baseline": {"lambda_gc": _mean([r["baseline"] for r in rows], "lambda_gc"),
                            "power": _mean([r["baseline"] for r in rows], "power"),
                            "n_causal_recovered": _mean(
                                [r["baseline"] for r in rows], "n_causal_recovered")},
               "rounds": {}, "qc": {
                   "histograms_identical": all(r["qc"]["histograms_identical"] for r in rows),
                   "max_abs_variance_change": max(r["qc"]["max_abs_variance_change"] for r in rows),
                   "max_abs_row_heterozygosity_change": max(
                       r["qc"]["max_abs_row_heterozygosity_change"] for r in rows)}}
        for strat in ("passive", "perturbed"):
            allr = [rr for r in rows for rr in r["rounds"][strat]]
            blk["rounds"][strat] = {
                "lambda_gc_mean": _mean(allr, "lambda_gc"),
                "power_mean": _mean(allr, "power"),
                "n_causal_recovered_mean": _mean(allr, "n_causal_recovered"),
                "n_causal_mean": _mean(allr, "n_causal"),
                "n_causal_in_qc_mask_mean": _mean(allr, "n_causal_in_qc_mask"),
                "n_off_causal_discoveries_mean": _mean(allr, "n_off_causal_discoveries"),
                "n_test_off_causal_mean": _mean(allr, "n_test_off_causal"),
                "off_causal_fraction_mean": _mean(allr, "off_causal_fraction"),
                "effect_agreement_mean": _mean(allr, "effect_agreement"),
                "min_p_causal_mean": _mean(allr, "min_p_causal")}
        out["cohorts"][name] = blk
        print("### %s  baseline lam %.3f power %.3f" % (
            name, blk["baseline"]["lambda_gc"], blk["baseline"]["power"]), flush=True)
        for strat in ("passive", "perturbed"):
            b = blk["rounds"][strat]
            print("   %-9s lam %.3f  power %.3f  offcausal %.2e  agree %.3f  minp %.1e" % (
                strat, b["lambda_gc_mean"], b["power_mean"], b["off_causal_fraction_mean"],
                b["effect_agreement_mean"], b["min_p_causal_mean"]), flush=True)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s71_protection_cost_chr{a.chrom}.json"), "w") as f:
        json.dump(payload, f, indent=2)
    print("done", payload["elapsed_s"], "s", flush=True)
