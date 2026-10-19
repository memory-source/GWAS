"""Common library for audit-governed multi-site GWAS experiments.

CPU-only, deterministic. Contains:
  - structured multi-site genotype simulator (parametric, Fst-based) with
    optional msprime-backed mode
  - standard GWAS utility pipeline (QC, PCA, marginal association, lambda_GC)
  - release-object construction (F1-F5 final tiers, M0-M3 intermediate tiers)
  - two membership-inference attack families + ancestry-only control attack
  - metrics (AUROC, TPR@FPR)
"""
from __future__ import annotations

import numpy as np


# ----------------------------------------------------------------------------
# Simulator
# ----------------------------------------------------------------------------
def _balding_nichols(base_freq, fst, rng, n_groups):
    """Draw ancestry-specific allele frequencies under the Balding-Nichols model."""
    base = np.clip(base_freq, 1e-3, 1 - 1e-3)
    if fst <= 0:
        return np.tile(base[None, :], (n_groups, 1))
    alpha = base * (1 - fst) / fst
    beta = (1 - base) * (1 - fst) / fst
    return rng.beta(alpha[None, :], beta[None, :], size=(n_groups, base.size))


def simulate_population(
    n_groups=3,
    n_per_group=700,
    n_snps=4000,
    fst=0.08,
    n_causal=10,
    causal_effect=0.35,
    stratification=0.6,
    ld_blocks=100,
    ld_size=6,
    seed=0,
):
    """Simulate ONE structured population that is then split into study members
    and non-members sharing the same ancestry and allele-frequency base.

    This is the controlled design the frontier requires: members and
    non-members differ only by study inclusion, not by population.
    """
    rng = np.random.default_rng(seed)
    base_freq = rng.uniform(0.05, 0.95, size=n_snps)
    group_freq = _balding_nichols(base_freq, fst, rng, n_groups)  # (K, m)

    ancestry = np.repeat(np.arange(n_groups), n_per_group)
    n = ancestry.size
    G = np.zeros((n, n_snps), dtype=np.int8)
    for k in range(n_groups):
        idx = np.where(ancestry == k)[0]
        G[np.ix_(idx, np.arange(n_snps))] = rng.binomial(
            2, group_freq[k][None, :], size=(idx.size, n_snps))

    if ld_blocks > 0 and ld_size > 1:
        n_blocks = min(ld_blocks, n_snps // ld_size)
        for b in range(n_blocks):
            cols = np.arange(b * ld_size, (b + 1) * ld_size)
            latent = rng.normal(0, 1, size=n)
            hi = latent > 0.5
            for c in cols:
                flip = rng.random(n) < 0.15
                G[hi & flip, c] = np.minimum(2, G[hi & flip, c] + 1)

    causal = rng.choice(n_snps, size=n_causal, replace=False)
    risk = G[:, causal].astype(np.float64) @ rng.normal(causal_effect, 0.05, size=n_causal)
    group_shift = (ancestry.astype(np.float64) - (n_groups - 1) / 2.0) * stratification
    logit = risk + group_shift - np.mean(risk + group_shift)
    prob = 1 / (1 + np.exp(-logit))
    case = (rng.random(n) < prob).astype(np.int8)

    return dict(G=G, ancestry=ancestry, group_freq=group_freq, base_freq=base_freq,
                case=case, causal=causal, n_groups=n_groups, n_snps=n_snps, seed=seed)


def split_population(pop, n_site_per_group=300, n_ref_per_group=300):
    """Split a simulated population into (study members, disjoint reference)
    per ancestry, so the two sets are ancestry-matched and population-matched.

    Returns two panel dicts. ``n_ref_per_group`` individuals drawn disjointly
    from the same population serve as the attacker's disjoint reference panel.
    """
    site, ref, _ = split_population3(pop, n_site_per_group, n_ref_per_group, 0)
    return site, ref


def split_population3(pop, n_site_per_group=300, n_ref_per_group=300,
                      n_nonmem_per_group=300):
    """Three-way disjoint split: study members, attacker reference panel, and
    an INDEPENDENT non-member target set (disjoint from BOTH the study and the
    attacker's reference).

    The non-member target panel is what the reviewer's disjoint-negative fix
    requires: scoring the attacker's own reference individuals as non-members
    is self-referential (they are pulled toward ref_af) and inflates the
    measured separation. ``n_nonmem_per_group=0`` reproduces the old (biased)
    behaviour for comparison.
    """
    anc = pop["ancestry"]
    G = pop["G"]
    case = pop["case"]

    def take(start, count):
        idx = []
        for k in range(pop["n_groups"]):
            cand = np.where(anc == k)[0]
            idx.append(cand[start:start + count])
        return np.concatenate(idx)

    site_idx = take(0, n_site_per_group)
    ref_idx = take(n_site_per_group, n_ref_per_group)
    nonmem_idx = (take(n_site_per_group + n_ref_per_group, n_nonmem_per_group)
                  if n_nonmem_per_group > 0 else np.array([], dtype=int))

    def sub(idx):
        return dict(G=G[idx], ancestry=anc[idx], case=case[idx],
                    group_freq=pop["group_freq"], base_freq=pop["base_freq"],
                    causal=pop["causal"], n_groups=pop["n_groups"],
                    n_snps=pop["n_snps"], seed=pop["seed"])

    return sub(site_idx), sub(ref_idx), sub(nonmem_idx)


def simulate_panel(
    n_groups=3,
    n_per_group=400,
    n_snps=4000,
    fst=0.08,
    n_causal=10,
    causal_effect=0.35,
    stratification=0.6,
    ld_blocks=100,
    ld_size=6,
    seed=0,
):
    """Backward-compatible entry point: returns a single combined panel."""
    return simulate_population(
        n_groups=n_groups, n_per_group=n_per_group, n_snps=n_snps, fst=fst,
        n_causal=n_causal, causal_effect=causal_effect,
        stratification=stratification, ld_blocks=ld_blocks, ld_size=ld_size,
        seed=seed,
    )


# ----------------------------------------------------------------------------
# Utility pipeline
# ----------------------------------------------------------------------------
def maf(G):
    """Minor allele frequency per SNP. NOTE: G.mean(0) is the ALT allele
    frequency times 2; a previous version treated it as the frequency directly
    and silently dropped every SNP with f>0.475."""
    f = G.mean(0) / 2.0
    return np.minimum(f, 1 - f)


def allele_freq(G):
    return G.mean(0) / 2.0


###--- QC mask ---###
def qc_mask(G, maf_floor=0.05, missing_rate=0.0, rng=None):
    """Return a boolean SNP mask. Genotypes here are complete; missing_rate
    simulates a QC summary statistic rather than masking the matrix."""
    m = maf(G)
    return (m >= maf_floor) & (m <= 1 - maf_floor)


def pca(G, k=10):
    """Top-k PCA on centered genotype matrix (returns eigenvalues, scores)."""
    X = G.astype(np.float64)
    X = X - X.mean(0, keepdims=True)
    # thin SVD via covariance in individual space is fine for n << m
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    eig = (S ** 2) / max(1, (X.shape[0] - 1))
    scores = U[:, :k] * S[:k]
    return eig[:k], scores, Vt[:k]


def genomic_control_lambda(pvals):
    chi2 = None
    from scipy.stats import chi2 as _chi2
    chi2 = _chi2.isf(np.clip(pvals, 1e-300, 1.0), 1)
    return np.median(chi2) / _chi2.ppf(0.5, 1)


def marginal_association(G, case, n_pcs=10, correct_pcs=True):
    """Chi-square association of each SNP with case/control, optionally
    correcting for the top principal components (population stratification)."""
    from scipy.stats import chi2 as _chi2

    n = G.shape[0]
    mask = qc_mask(G)
    X = G[:, mask].astype(np.float64)
    y = case.astype(np.float64)

    if correct_pcs and n_pcs > 0:
        _, scores, _ = pca(G, k=n_pcs)
        design = np.column_stack([np.ones(n), scores])
        # residualize both X and y on the PCs
        beta_x, *_ = np.linalg.lstsq(design, X, rcond=None)
        Xr = X - design @ beta_x
        beta_y, *_ = np.linalg.lstsq(design, y, rcond=None)
        yr = y - design @ beta_y
    else:
        Xr, yr = X - X.mean(0), y - y.mean()

    # per-SNP chi-square (1 df) for mean difference between groups
    xc = Xr - Xr.mean(0)
    yc = yr - yr.mean()
    denom = np.sqrt((xc ** 2).sum(0) * (yc ** 2).sum())
    denom[denom == 0] = np.inf
    r = (xc * yc[:, None]).sum(0) / denom
    r = np.clip(r, -1 + 1e-12, 1 - 1e-12)
    stat = r ** 2 * n
    pvals = _chi2.sf(stat, 1)
    return mask, pvals


# ----------------------------------------------------------------------------
# Release objects
# ----------------------------------------------------------------------------
def release_objects(panel, n_pcs=10, rng_seed=0, study_only=True):
    """Compute the quantities a real multi-site GWAS workflow publishes.

    ``study_only=True`` estimates PCA/allele frequencies from the study panel
    only (which is what a study actually releases); the attacker's reference
    is supplied separately.
    """
    G = panel["G"]
    case = panel["case"]
    mask, pvals = marginal_association(G, case, n_pcs=n_pcs, correct_pcs=True)
    kept_idx = np.where(mask)[0]
    p_full = np.full(G.shape[1], np.nan)
    p_full[kept_idx] = pvals
    af = allele_freq(G)
    eig, scores, loadings = pca(G, k=n_pcs)
    lam = genomic_control_lambda(pvals)
    rng = np.random.default_rng(rng_seed)
    missing = rng.uniform(0.0, 0.02, size=G.shape[1])
    # per-site (per-stratum) allele frequencies: the realistic multi-site
    # release unit of a federated meta-analysis.
    anc = panel.get("ancestry")
    af_by_group = {}
    if anc is not None:
        for k in np.unique(anc):
            af_by_group[int(k)] = allele_freq(G[anc == k])
    return dict(
        mask=mask, kept_idx=kept_idx, pvals=p_full,
        af=af, af_by_group=af_by_group, eig=eig,
        loadings=loadings, scores=scores,
        lambda_gc=lam, missing=missing,
    )


# ----------------------------------------------------------------------------
# Attack families
# ----------------------------------------------------------------------------
def _hwe_geno_loglik(g, f):
    """log P(genotype g | allele freq f) under HWE; g in {0,1,2}."""
    f = np.clip(f, 1e-6, 1 - 1e-6)
    p0 = (1 - f) ** 2
    p1 = 2 * f * (1 - f)
    p2 = f ** 2
    out = np.where(g == 0, np.log(p0), np.where(g == 1, np.log(p1), np.log(p2)))
    return out


def af_lrt_attack(target_G, study_af, ref_af, mask=None):
    """Homer/Sankararaman-style likelihood-ratio membership score.

    Score = sum_j [ log P(g_ij | study_af_j) - log P(g_ij | ref_af_j) ].
    Members (whose genotype contributed to study_af) score higher.
    target_G: (n_targets, m) or (m,). study_af/ref_af: (m,).
    """
    g = np.atleast_2d(target_G)
    if mask is not None:
        g = g[:, mask]
        study_af = study_af[mask]
        ref_af = ref_af[mask]
    ll_study = _hwe_geno_loglik(g, study_af[None, :])   # (n, m)
    ll_ref = _hwe_geno_loglik(g, ref_af[None, :])
    return (ll_study - ll_ref).sum(1)


def beacon_lrt_attack(target_G, study_G, mask=None, max_snps=2000, rng=None):
    """Shringarpure/Bustamante-style beacon likelihood-ratio attack.

    Two parameters:
      N  = beacon size = number of individuals in the study (determines the
           false-positive rate p0 of a "no" answer under the null);
      f  = *background* (population/reference) alt-allele frequency, NOT the
           study's own - the attacker estimates background from an external
           reference panel.
    The attacker queries only positions where the target carries >=1 alt allele
    (the informative set).
    """
    rng = np.random.default_rng(0 if rng is None else rng)
    g = np.atleast_2d(target_G)
    study_G = np.atleast_2d(study_G)
    if mask is not None:
        g = g[:, mask]
        study_G = study_G[:, mask]
    n = study_G.shape[0]
    resp = (study_G > 0).any(0).astype(np.float64)
    f = float(np.clip(np.mean(study_G) / 2.0, 1e-3, 1 - 1e-3))
    p0 = float(np.clip(1 - (1 - f) ** (2 * n), 1e-6, 1 - 1e-6))

    carries = g > 0
    if max_snps is not None and carries.shape[1] > max_snps:
        sel = rng.choice(carries.shape[1], size=max_snps, replace=False)
        carries = carries[:, sel]
        resp = resp[sel]

    yes = resp[None, :]
    ll1 = yes * np.log(1.0) + (1 - yes) * np.log(1e-12)
    ll0 = yes * np.log(p0) + (1 - yes) * np.log(1 - p0)
    contrib = (ll1 - ll0) * carries
    return contrib.sum(1)


def anc_only_attack(target_G, study_loadings, study_mean, mask):
    """Confound upper bound: score a target by its projection onto the
    STUDY-fitted PC basis, using no allele-frequency information.

    The target is centered by the (public) study mean, projected onto the
    study's principal axes, and scored by proximity to the study centroid
    (which is the origin in the study-centered basis). This is the
    ancestry/structure channel an attacker could exploit without any released
    AF; it is the correct control (an earlier helper was a no-op: it used each
    panel's own -- centered -- PCA scores, so the centroid was 0 and the study
    members were never scored).
    """
    X = target_G[:, mask].astype(np.float64) - study_mean[None, :]
    proj = (X @ study_loadings[:, mask].T) / np.sqrt(max(1, mask.sum()))
    return -np.linalg.norm(proj, axis=1)


def ancestry_only_attack(target_scores, ref_scores):
    """Deprecated no-op control retained for callers that pass PCA scores.
    Use ``anc_only_attack`` with study-fitted loadings instead."""
    mu = ref_scores.mean(0)
    return -np.linalg.norm(target_scores - mu[None, :], axis=1)


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------
def auroc(scores, labels):
    """Rank-sum AUROC with midrank tie handling. For all-constant scores this
    returns 0.5 (not 0.0), so a degenerate indicator column does not masquerade
    as a far-below-chance classifier."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    s = np.concatenate([pos, neg])
    order = np.argsort(s, kind="mergesort")
    sorted_s = s[order]
    _, inv, counts = np.unique(sorted_s, return_inverse=True, return_counts=True)
    csum = np.cumsum(counts)
    avg = (csum - counts + 1 + csum) / 2.0          # average rank within each tie group
    ranks = np.empty(len(s), dtype=float)
    ranks[order] = avg[inv]
    r_pos = ranks[: pos.size].sum()
    return (r_pos - pos.size * (pos.size + 1) / 2) / (pos.size * neg.size)


def tpr_at_fpr(scores, labels, fpr_target):
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    thresholds = np.sort(np.unique(scores))[::-1]
    best = 0.0
    for t in thresholds:
        pred = scores >= t
        tp = np.sum(pred & (labels == 1))
        fp = np.sum(pred & (labels == 0))
        tpr = tp / max(1, labels.sum())
        fpr = fp / max(1, (labels == 0).sum())
        if fpr <= fpr_target:
            best = max(best, tpr)
    return best


def stratified_permutation_auroc(scores, labels, strata, n_perm=200, seed=0):
    """AUC under within-stratum label permutation: the individual-membership
    component after removing ancestry/cohort confound."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_perm):
        perm = labels.copy()
        for s in np.unique(strata):
            idx = np.where(strata == s)[0]
            perm[idx] = rng.permutation(perm[idx])
        out.append(auroc(scores, perm))
    return float(np.mean(out)), float(np.std(out))
