"""S61 - selection-input spectral tomography (deep-search direction).

The passive paper measures leakage from a FIXED release. This script measures the
ACTIVE question of `research-plan/20261016_deep_search/DEEP_DIRECTIONS_ZH.md`:
a cooperating analysis partner may vary its own AF-preserving input across re-fits
of a joint PCA and reads only the final truncated output each time.

Model (fixed public center and scale, matching the direction document):
  honest scatter   C_H = Z_H^T Z_H          (Z_H = (honest - mu)/scale)
  partner scatter  K_t = Z_P(t)^T Z_P(t)    (Z_P(t) = (own_t - mu)/scale, known)
  joint eigh       (C_H + K_t) V_t = V_t Lambda_t
  final output     (V_t, Lambda_t)  =>  the partner forms
                     Y_t = V_t Lambda_t - K_t V_t  =  C_H V_t     (exact identity)
  stacking q       Q = [V_1..V_q], Y = [Y_1..Y_q],  C_H Q = Y
                   =>  recovered measure  C_H P_Q, P_Q the projector on span(Q).

Input change (AF-preserving): exchange one contiguous SNP block between two of the
partner's own rows. Every per-SNP genotype histogram, AF, variance, HWE count and
per-SNP missing count is unchanged; only the individual correspondence changes.

Roles are drawn from real genotypes WITHOUT phenotype ascertainment; a phenotype
exists only for the GWAS-power row (a simulated case/control with known causal
markers), and that row fits the association on the honest study alone.

Target directions: the passive first release defines a span P0. The attacker's
independent shadow panel supplies candidate directions; x_in is the shadow
direction with the SMALLEST residual after projecting out span(P0) (already
reached by the fixed release) and x_out the LARGEST (outside the fixed release).
Both are defined by construction against the passive span, not asserted.

Strategies (x-axis = number of final outputs q in {1,2,4,8}):
  passive  -- the partner input never changes (all q outputs identical, span = k)
  random   -- one random AF-preserving block swap per release
  adaptive -- the swap is chosen to align the *predicted* joint top-k with x_out and
              add direction not in span(Q_prior); candidate ranking uses ONLY the
              partner matrix, the independent shadow panel, past released
              eigenvectors and the public center/scale (never the honest private
              genotype matrix).

Output: results/raw/s61_spectral_tomography.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s48_real_structure as s48
import s51_real_lowstruct as s51
from common import (allele_freq, af_lrt_attack, auroc, tpr_at_fpr,
                    genomic_control_lambda, marginal_association)
from real_multi import SUPERPOP

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "raw")
QS = (1, 2, 4, 8)
STRATEGIES = ("passive", "random", "adaptive")
K_TOP = 10
BLOCK_LEN = 20
N_CANDIDATES = 24
N_PERM = 2000
N_BOOT = 10000
# |z| bound used only for the no-release controls: the EAS/3pop panels are
# MAF-filtered on the full 2504-sample panel and carry globally-common markers that
# are monomorphic inside the cohort (cohort-mu == 0), where the scale clip (1e-6)
# would give a random span a spurious anti-membership signal.
Z_CLIP = 5.0
# per-GROUP caps; the sum (420) must stay under the smallest group in any cohort
# (EUR has 503 chr22 rows at the 1.5 kb thinning). Roles are disjoint per group.
CAPS = {"honest": 120, "partner": 120, "nonmember": 80,
        "reference": 60, "shadow": 40}
# per-GROUP caps for the honest-sample-size ladder (largest is used for the
# sample-size arm of the sensitivity check; EUR chr22 has 503 rows at 1.5 kb).
CAPS_BIG = {"honest": 180, "partner": 180, "nonmember": 80,
            "reference": 60, "shadow": 40}


# ----------------------------------------------------------------- helpers ---
def split_roles(groups, caps, seed):
    """Disjoint roles per group label, so every role is structure-matched. Falls
    over to an error if the pool cannot supply the requested counts (never a silent
    truncation, which would change a role's size without saying so)."""
    rng = np.random.default_rng(seed)
    roles = {r: [] for r in caps}
    for g in np.unique(groups):
        idx = rng.permutation(np.where(groups == g)[0])
        start = 0
        for r, cap in caps.items():
            if start + cap > idx.size:
                raise ValueError(
                    f"group {g} has {idx.size} rows, need {start + cap} for the "
                    f"requested role caps {caps}")
            roles[r].append(idx[start:start + cap])
            start += cap
    return {r: np.concatenate(v) for r, v in roles.items()}


def zscore(rows, mu, scale):
    return (rows.astype(np.float64) - mu[None, :]) / scale[None, :]


def joint_baseline(honest, partner):
    """Public center/scale from the baseline joint (honest + partner-initial)."""
    mu = np.concatenate([honest, partner]).mean(axis=0)
    scale = np.sqrt(np.clip(mu * (1 - mu / 2.0), 1e-6, None))
    return mu, scale


def topk_gram(Z, k):
    """Top-k eigenvectors of Z^T Z for one stacked centered matrix."""
    Gm = Z @ Z.T
    ev, U = np.linalg.eigh(Gm)
    idx = np.argsort(ev)[::-1][:k]
    lam = ev[idx]
    S = np.sqrt(np.clip(lam, 1e-12, None))
    V = (Z.T @ U[:, idx]) / S[None, :]
    return lam, V


def topk(Zs, k):
    """Top-k eigenvectors/values of sum_i Z_i^T Z_i (several centered blocks)."""
    if len(Zs) == 1:
        return topk_gram(Zs[0], k)
    return topk_gram(np.vstack(Zs), k)


def swap_block(own, pair, block, block_len=BLOCK_LEN):
    """AF-preserving exchange of one contiguous SNP block between two own rows.
    Returns (new_matrix, n_genotypes_actually_changed) where the count is the
    number of MATRIX ENTRIES that change (two per differing marker position)."""
    out = own.copy()
    c0 = block * block_len
    c1 = min(c0 + block_len, own.shape[1])
    diff_positions = int((own[pair[0], c0:c1] != own[pair[1], c0:c1]).sum())
    out[np.ix_(pair, np.arange(c0, c1))] = own[np.ix_(pair[::-1], np.arange(c0, c1))]
    return out, 2 * diff_positions


def hist_sig(own):
    return np.sort(own, axis=0)


def sim_phenotype(site, seed, n_causal=8):
    """Mechanism-clear simulated case/control on the study rows: a few causal
    markers with additive effects. Used for the GWAS-power row only."""
    rng = np.random.default_rng(seed + 900)
    m = site.shape[1]
    causal = rng.choice(m, size=n_causal, replace=False)
    eff = rng.normal(0.0, 0.35, size=n_causal)
    risk = site[:, causal].astype(float) @ eff
    risk = risk - risk.mean()
    prob = 1.0 / (1.0 + np.exp(-risk))
    return (rng.random(site.shape[0]) < prob).astype(np.int8), causal


def perm_pvalue(ev_scores, y_ev, strata_ev, n_perm, seed):
    """P(permuted AUROC >= observed), permuting labels within strata on the
    evaluated rows only (scores fixed)."""
    obs = auroc(ev_scores, y_ev)
    rng = np.random.default_rng(seed)
    groups = [np.where(strata_ev == g)[0] for g in np.unique(strata_ev)]
    ge = 0
    for _ in range(n_perm):
        perm = y_ev.copy()
        for idx in groups:
            perm[idx] = rng.permutation(perm[idx])
        ge += int(auroc(ev_scores, perm) >= obs)
    return (ge + 1) / (n_perm + 1)


# ----------------------------------------------------------------- blocks ---
def witness(G, groups, seed=0, caps=None):
    """B0: identity residual Y_t = C_H V_t, AF-preservation audit, and the
    individual-level movement, all measured over the SAME 12 releases."""
    roles = split_roles(groups, caps or CAPS, seed)
    honest, partner = G[roles["honest"]], G[roles["partner"]]
    mu, scale = joint_baseline(honest, partner)
    Zh = zscore(honest, mu, scale)
    C_H = Zh.T @ Zh
    sig0 = hist_sig(partner)
    af0, var0 = allele_freq(partner), partner.astype(float).var(0)
    het0 = (partner == 1).mean(1)
    n_blocks = G.shape[1] // BLOCK_LEN
    rng = np.random.default_rng(seed + 7)
    max_res = max_af = max_var = max_het = 0.0
    hist_ok = True
    changed_counts = []
    for _ in range(12):
        pair = rng.choice(partner.shape[0], 2, replace=False)
        own, nchg = swap_block(partner, pair, int(rng.integers(0, n_blocks)))
        changed_counts.append(nchg)
        hist_ok = hist_ok and bool(np.array_equal(hist_sig(own), sig0))
        max_af = max(max_af, float(np.max(np.abs(allele_freq(own) - af0))))
        max_var = max(max_var, float(np.max(np.abs(own.astype(float).var(0) - var0))))
        max_het = max(max_het, float(np.max(np.abs((own == 1).mean(1) - het0))))
        Zp = zscore(own, mu, scale)
        K = Zp.T @ Zp
        lam, V = topk([Zh, Zp], K_TOP)
        Y = V * lam[None, :] - K @ V
        max_res = max(max_res, float(np.max(np.abs(C_H @ V - Y))))
    return {
        "snps": int(G.shape[1]), "honest_rows": int(len(roles["honest"])),
        "partner_rows": int(len(roles["partner"])), "releases": 12,
        "block_len": BLOCK_LEN,
        "selected_rows_per_release": 2,
        "genotypes_changed_per_release": changed_counts,
        "changed_fraction_of_partner_rows": 2.0 / len(roles["partner"]),
        "max_eigen_equation_residual": max_res,
        "histograms_identical_all_releases": hist_ok,
        "max_abs_af_change": max_af, "max_abs_var_change": max_var,
        "max_abs_row_heterozygosity_change": max_het,
    }


def adaptive_input(partner, x_target, Q_prior, mu, scale, Zsh, n_blocks, rng, k=K_TOP):
    """Greedy block swap: pick the candidate whose predicted joint top-k both aligns
    with x_target and adds direction not already in span(Q_prior). Uses partner data,
    the independent shadow panel, past released eigenvectors and the public
    center/scale; never the honest private genotype matrix."""
    if Q_prior:
        Qp = np.concatenate(Q_prior, axis=1)
        # same relative rank cutoff as the evaluation, so the prior span used for
        # the gain term matches the span the released outputs actually define
        Uq, sq, Vtq = np.linalg.svd(Qp, full_matrices=False)
        keepp = sq > 1e-9 * (sq[0] if sq.size else 1.0)
        Bp = Uq[:, keepp]
    else:
        Bp = None
    xt = x_target / max(1e-12, np.linalg.norm(x_target))

    def score(Vc):
        align = float(np.linalg.norm(Vc.T @ xt))
        if Bp is None:
            return align
        covered = Bp @ (Bp.T @ Vc)                     # projection onto span(Q_prior)
        gain = float(np.linalg.norm(Vc - covered))
        return align * gain

    best, best_score = partner, -1.0
    for _ in range(N_CANDIDATES):
        pair = rng.choice(partner.shape[0], 2, replace=False)
        cand, _ = swap_block(partner, pair, int(rng.integers(0, n_blocks)))
        Zc = zscore(cand, mu, scale)
        _, Vc = topk([Zsh, Zc], k)
        sc = score(Vc)
        if sc > best_score:
            best_score, best = sc, cand
    return best


def _target_dirs(Zsh, P0):
    """x_in / x_out: the shadow directions with the smallest / largest residual
    after projecting out the passive first-release span P0 (orthonormal)."""
    _, Vsh = topk([Zsh], min(K_TOP + 1, Zsh.shape[0] - 1))
    res = []
    for j in range(Vsh.shape[1]):
        x = Vsh[:, j]
        resid = np.linalg.norm(x - P0 @ (P0.T @ x))
        res.append((resid, x))
    res.sort(key=lambda t: t[0])
    return res[0][1], res[-1][1], float(res[0][0]), float(res[-1][0])


def metrics(Q, Y, x_in, x_out, honest, nonmem, groups, mu, scale, y, cal, seed,
            dedup=True, do_perm=True, C_H=None):
    """Coverage and membership metrics for one stacked output (Q, Y).

    `dedup=True` applies the relative rank cutoff to the visible basis and the
    pseudo-inverse (used for the stacked comparisons, where passive repeats must
    not inflate the span). `dedup=False` uses the full pseudo-inverse, used only for
    the single-release rows where Q is expected to be near-square and full rank.
    `C_H` is only needed for the coverage denominators (evaluator-side); the attack
    metrics do not read it."""
    Uq, sq, Vtq = np.linalg.svd(Q, full_matrices=False)
    if dedup:
        keep = sq > 1e-9 * (sq[0] if sq.size else 1.0)
        Bq = Uq[:, keep]
        sqk = sq[keep]
        Qpinv = (Vtq[keep].T / sqk[None, :]) @ Uq[:, keep].T
    else:
        keep = sq > 0
        Bq = Uq
        Qpinv = np.linalg.pinv(Q)
    norm_CH = float(np.linalg.norm(C_H)) if C_H is not None else 1.0
    visible = float(np.linalg.norm(C_H @ Bq) / norm_CH) if C_H is not None else float("nan")
    visible_energy = visible ** 2
    scatter_err = (float(np.linalg.norm(Y @ Qpinv - C_H) / norm_CH)
                   if C_H is not None else float("nan"))

    def tgt(x):
        if C_H is None:
            return float("nan")
        true = float(x @ C_H @ x)
        est = float(x @ Y @ (Qpinv @ x))
        return abs(est - true) / max(1e-12, abs(true))

    def tres(x):
        return float(np.linalg.norm(x - Q @ (Qpinv @ x)) / max(1e-12, np.linalg.norm(x)))

    # membership from the range of the released outputs Y (M x q); cheap SVD, and
    # passive repeats naturally collapse under the same relative rank cutoff.
    Wy, sy, _ = np.linalg.svd(Y, full_matrices=False)
    keepY = sy > 1e-9 * (sy[0] if sy.size else 1.0)
    B = Wy[:, keepY]
    rank_recovered = int(keepY.sum())
    X = zscore(np.vstack([honest, nonmem]), mu, scale)
    s = np.linalg.norm(X @ B, axis=1)
    h = np.asarray(cal, bool)
    sgn = 1.0 if auroc(s[h], y[h]) >= 0.5 else -1.0
    s_ev, y_ev = sgn * s[~h], y[~h]
    rec_auc = float(auroc(s_ev, y_ev))
    tpr01 = float(tpr_at_fpr(s_ev, y_ev, 0.01))
    tpr05 = float(tpr_at_fpr(s_ev, y_ev, 0.05))
    p_perm = (perm_pvalue(s_ev, y_ev, groups[~h], N_PERM, seed + 3)
              if do_perm else float("nan"))
    return {"rank_Q": int(np.linalg.matrix_rank(Q)),
            "rank_visible": int(keep.sum()),
            "rank_recovered": rank_recovered,
            "min_sv": float(sq[-1]) if sq.size else 0.0,
            "visible_share": visible,
            "visible_energy": visible_energy,
            "scatter_recon_err": scatter_err,
            "target_resid_inside": tres(x_in),
            "target_resid_outside": tres(x_out),
            "target_err_inside": tgt(x_in),
            "target_err_outside": tgt(x_out),
            "membership_auc": rec_auc,
            "membership_tpr_fpr01": tpr01,
            "membership_tpr_fpr05": tpr05,
            "membership_perm_p": p_perm}


def _membership_auc_from_operator(R, honest, nonmem, mu, scale, y, cal, sign=None,
                                  rel_cutoff=1e-9, z_clip=None):
    """Membership AUROC (and TPR@5%FPR) from an arbitrary recovered operator R
    (range = span). Pass `sign` to reuse an orientation fixed elsewhere (needed when
    conditioning removes the calibration members and leaves that half uninformative).
    `z_clip` bounds |z| for the rows entering the score; it is used only for the
    no-release controls, where the EAS/3pop panels carry globally-common markers that
    are monomorphic inside the cohort (mu == 0) and would give a random span a
    spurious anti-membership signal through the scale clip. Returns (auc, tpr05, sign)."""
    U, s, _ = np.linalg.svd(R, full_matrices=False)
    keep = s > rel_cutoff * (s[0] if s.size else 1.0)
    B = U[:, keep]
    X = zscore(np.vstack([honest, nonmem]), mu, scale)
    if z_clip is not None:
        X = np.clip(X, -z_clip, z_clip)
    d = np.linalg.norm(X @ B, axis=1)
    h = np.asarray(cal, bool)
    sgn = sign if sign is not None else (1.0 if auroc(d[h], y[h]) >= 0.5 else -1.0)
    s_ev, y_ev = sgn * d[~h], y[~h]
    return float(auroc(s_ev, y_ev)), float(tpr_at_fpr(s_ev, y_ev, 0.05)), sgn


def af_baselines(honest, nonmem, ref, partner, shadow, study_groups, y, cal, seed):
    """AF-LRT baselines under the SAME allowed information the spectral attack has.

    The attacker is granted the honest cohort's released AF, the partner matrix and
    the independent shadow panel (both must not contain the scored targets). Three
    baselines:
      af_ref       -- background frequency from the disjoint reference alone
      af_pooled    -- background frequency pooled over reference+partner+shadow
      af_cond_pool -- composition-balanced removal of the known calibration members
                      from the released AF, then the pooled background frequency
    Each held-out-half oriented (sign fixed on the unconditioned channel's
    calibration half)."""
    af_s = allele_freq(honest)
    bg_pooled = allele_freq(np.vstack([ref, partner, shadow]))
    h = np.asarray(cal, bool)
    outs = []
    for name, bg in (("ref", allele_freq(ref)), ("pooled", bg_pooled)):
        af_raw = np.concatenate([af_lrt_attack(honest, af_s, bg),
                                 af_lrt_attack(nonmem, af_s, bg)])
        outs.append(s48._orient_auroc(af_raw, y, np.random.default_rng(seed + 31),
                                      cal_mask=cal))
        if name == "pooled":
            sgn_af = 1.0 if auroc(af_raw[h], y[h]) >= 0.5 else -1.0
    af_auc = outs[0]
    af_pooled_auc = outs[1]
    # conditioned (composition-balanced removal) on the pooled background
    N = len(honest)
    gh = np.asarray(study_groups)
    h_study = h[:N]
    groups_present = np.unique(gh)
    quota = min(int((h_study & (gh == g)).sum()) for g in groups_present)
    bal = np.concatenate([np.where(h_study & (gh == g))[0][:quota] for g in groups_present])
    if quota > 0 and len(bal) < N:
        num2 = 2.0 * N * af_s - honest[bal].sum(0).astype(float)
        af_cond = np.clip(num2 / (2.0 * (N - len(bal))), 1e-6, 1 - 1e-6)
        cond_raw = np.concatenate([af_lrt_attack(honest, af_cond, bg_pooled),
                                   af_lrt_attack(nonmem, af_cond, bg_pooled)])
        af_cond_auc = float(auroc(sgn_af * cond_raw[~h], y[~h]))
    else:
        af_cond_auc = float("nan")
    return af_auc, af_pooled_auc, af_cond_auc, int(len(bal))


def cell(G, groups, seed, k=K_TOP, block_len=BLOCK_LEN, caps=None, roles=None, qs=QS):
    if roles is None:
        roles = split_roles(groups, caps or CAPS, seed)
    honest, partner = G[roles["honest"]], G[roles["partner"]]
    nonmem, ref, shadow = G[roles["nonmember"]], G[roles["reference"]], G[roles["shadow"]]
    mu, scale = joint_baseline(honest, partner)
    Zh = zscore(honest, mu, scale)
    Zsh = zscore(shadow, mu, scale)
    C_H = Zh.T @ Zh
    n_blocks = G.shape[1] // block_len

    y = np.concatenate([np.ones(len(honest)), np.zeros(len(nonmem))])
    cal = np.random.default_rng(seed + 31).random(len(y)) < 0.5
    grp_rows = np.concatenate([np.asarray(groups)[roles["honest"]],
                               np.asarray(groups)[roles["nonmember"]]])
    honest_groups = np.asarray(groups)[roles["honest"]]
    af_auc, af_pooled_auc, af_cond_auc, n_bal = af_baselines(
        honest, nonmem, ref, partner, shadow, honest_groups, y, cal, seed)

    # The passive first release is a PUBLIC baseline every strategy shares: the
    # attacker forms one fixed-input release, defines the target direction against
    # its span, and that release also counts as output 1. Each strategy then emits
    # further outputs always from a fresh partner input, so output q costs q
    # releases for every strategy (adaptive replaces its input, it does not spend
    # a free extra one).
    Zp0 = zscore(partner, mu, scale)
    lam0, V0 = topk([Zh, Zp0], k)
    P0 = np.linalg.qr(V0)[0]
    x_in, x_out, resid_in, resid_out = _target_dirs(Zsh, P0)

    # conditioned single-release spectral baseline: subtract the known calibration
    # members' scatter contribution from the release, then score its range. Uses
    # only information the plan grants the attacker. The orientation is fixed on the
    # UNCONDITIONED single release (the conditioning removes the calibration members
    # and leaves that half uninformative), exactly as the conditioned-AF baseline.
    h = np.asarray(cal, bool)
    K0m = Zp0.T @ Zp0
    Y0 = V0 * lam0[None, :] - K0m @ V0
    _, _, sgn_spec = _membership_auc_from_operator(Y0, honest, nonmem, mu, scale, y, cal)
    Zk = Zh[h[:len(honest)]]                    # calibration honest members' z
    Y0_cond = Y0 - Zk.T @ (Zk @ V0)
    spec_cond_auc, spec_cond_tpr, _ = _membership_auc_from_operator(
        Y0_cond, honest, nonmem, mu, scale, y, cal, sign=sgn_spec)

    # no-release control: the attacker holds the KNOWN calibration members' full
    # genotypes but has NO release. Score projection onto (i) their own top-k PCA
    # span, (ii) their top-(qk) span, (iii) a random d-dimensional span. If these
    # are at or below chance, the release-derived span is what carries the signal.
    # The rows are |z|-clipped: the EAS/3pop panels are MAF-filtered on the full
    # 2504-sample panel, so they carry globally-common markers that are monomorphic
    # inside the cohort (cohort-mu == 0); the scale clip would give a random span a
    # spurious anti-membership signal from those markers (see the results md).
    no_release = {}
    for tag, d in (("cal_pca_k", k), ("cal_pca_qk", min(max(qs) * k, Zk.shape[0] - 1)),
                   ("random_k", k)):
        if tag == "random_k":
            R = np.random.default_rng(seed + 71).standard_normal((Zk.shape[1], d))
        else:
            _, Vk = topk([Zk], d)
            R = Vk
        auc_c, _, _ = _membership_auc_from_operator(R, honest, nonmem, mu, scale, y, cal,
                                                    z_clip=Z_CLIP)
        no_release[tag] = auc_c

    out = {"af_auc": af_auc, "af_pooled_auc": af_pooled_auc,
           "af_conditioned_auc": af_cond_auc, "n_conditioned_members": n_bal,
           "spec_conditioned_auc": spec_cond_auc, "spec_conditioned_tpr05": spec_cond_tpr,
           "no_release_controls": no_release,
           "target_resid_inside": resid_in, "target_resid_outside": resid_out,
           "strategies": {}, "precision": {}}
    for strategy in STRATEGIES:
        rng = np.random.default_rng(seed + 101)
        # output 1 is the shared passive release (full pseudo-inverse: it is the
        # single near-square release, not a stacked comparison)
        Qs, Ys = [V0], [Y0]
        # per-release released quantities, kept for the export-precision probe
        chunks = [(V0.copy(), lam0.copy(), K0m.copy())]
        own_seq = [partner]                       # partner input used at each release
        per_q = {"1": metrics(V0, Y0, x_in, x_out, honest, nonmem,
                              grp_rows, mu, scale, y, cal, seed, dedup=False, C_H=C_H)}
        for t in range(1, max(qs)):
            if strategy == "passive":
                own = partner
            elif strategy == "random":
                own, _ = swap_block(partner, rng.choice(partner.shape[0], 2, replace=False),
                                    int(rng.integers(0, n_blocks)), block_len)
            else:
                own = adaptive_input(partner, x_out, Qs, mu, scale, Zsh, n_blocks, rng)
            Zp = zscore(own, mu, scale)
            K = Zp.T @ Zp
            lam, V = topk([Zh, Zp], k)
            Y = V * lam[None, :] - K @ V
            Qs.append(V); Ys.append(Y)
            chunks.append((V.copy(), lam.copy(), K.copy()))
            own_seq.append(own.copy())
            q = t + 1
            if q in qs:
                Q = np.concatenate(Qs, axis=1); Y = np.concatenate(Ys, axis=1)
                per_q[str(q)] = metrics(Q, Y, x_in, x_out, honest, nonmem,
                                        grp_rows, mu, scale, y, cal, seed, C_H=C_H)
        out["strategies"][strategy] = per_q
        if strategy == "random":
            if qs == QS:
                out["precision"] = _precision_probe(
                    C_H, chunks, (V0, lam0, K0m, Zk), af_cond_auc, x_in, x_out, honest,
                    nonmem, grp_rows, mu, scale, y, cal, seed)
                out["sci_cost"] = sci_cost(honest, own_seq, seed, mu, scale)
    return out


def _precision_probe(C_H, chunks, spec0, af_cond, x_in, x_out, honest, nonmem, groups,
                     mu, scale, y, cal, seed):
    """Re-form the published quantities at reduced numeric precision (float32, 6 and
    4 decimals) and re-run coverage, membership, and the paired gap against a
    PRECISION-MATCHED strong passive baseline (the stronger of the conditioned AF
    and the conditioned single release re-formed at the same precision). Only the
    RELEASED quantities (V_t, Lambda_t) are quantized; the partner's own scatter is
    known exactly."""
    V0, lam0, K0m, Zk = spec0
    h = np.asarray(cal, bool)

    def reform(V, lam, K):
        # Lambda and V are the RELEASED quantities (quantized); K is the partner's
        # own scatter, known exactly.
        return V * lam[None, :] - K @ V

    probe = {}
    for tag, transform in (("float64", lambda A: A),
                           ("float32", lambda A: A.astype(np.float32).astype(np.float64)),
                           ("round6", lambda A: np.round(A, 6)),
                           ("round4", lambda A: np.round(A, 4))):
        # precision-matched strong passive baseline (single release, conditioned)
        V0r, l0r = transform(V0), transform(lam0)
        Y0r = reform(V0r, l0r, K0m)
        _, _, sgn_p = _membership_auc_from_operator(Y0r, honest, nonmem, mu, scale, y, cal)
        Y0c = Y0r - Zk.T @ (Zk @ V0r)
        spec_p, _, _ = _membership_auc_from_operator(Y0c, honest, nonmem, mu, scale, y, cal,
                                                     sign=sgn_p)
        strong = max(af_cond, spec_p)
        Qs, Ys = [], []
        for V, lam, K in chunks:
            Vr, lr = transform(V), transform(lam)
            Ys.append(reform(Vr, lr, K))
            Qs.append(Vr)
        Q = np.concatenate(Qs, axis=1); Y = np.concatenate(Ys, axis=1)
        m = metrics(Q, Y, x_in, x_out, honest, nonmem, groups, mu, scale,
                    y, cal, seed, do_perm=False, C_H=C_H)
        probe[tag] = {"scatter_recon_err": m["scatter_recon_err"],
                      "min_sv": m["min_sv"], "rank_visible": m["rank_visible"],
                      "visible_share": m["visible_share"],
                      "membership_auc": m["membership_auc"],
                      "spec_cond_precision_auc": spec_p,
                      "gap_vs_strongest_passive": m["membership_auc"] - strong}
    return probe


def gwaspower(honest, seed):
    """Scientific-cost baseline: GWAS signal recovery on the HONEST STUDY ALONE
    with the honest study's OWN top-10 PCs as covariates."""
    case, causal = sim_phenotype(honest, seed)
    mask, p = marginal_association(honest, case, n_pcs=10, correct_pcs=True)
    hit = np.zeros(honest.shape[1], bool)
    hit[causal] = True
    hit &= mask
    hit_kept = hit[mask]
    n_hit = int(hit_kept.sum())
    return {"n_study": int(honest.shape[0]), "n_causal_in_mask": n_hit,
            "min_p_causal": float(np.min(p[hit_kept])) if n_hit else None,
            "lambda_gc": float(genomic_control_lambda(p))}


def sci_cost(honest, partner_seq, seed, zh_mean, zh_scale, k=K_TOP):
    """Per-release scientific cost: run the honest study's GWAS using the covariates
    from each round's JOINT PCA of (honest, partner_t), and compare lambda_gc and the
    causal-marker min-p against the honest-only baseline. `partner_seq` is the list
    of partner inputs the random strategy used at each release."""
    case, causal = sim_phenotype(honest, seed)
    Zh = zscore(honest, zh_mean, zh_scale)
    rows = []
    for own in partner_seq:
        Zp = zscore(own, zh_mean, zh_scale)
        lam, V = topk([Zh, Zp], k)          # V is (M, k)
        # per-individual joint PC scores on the honest rows (study-side covariates)
        scores = Zh @ V                                     # (n_honest, k)
        mask, p = marginal_association(honest, case, n_pcs=k, correct_pcs=False)
        # replace the auto-PC covariates with the joint-PCA scores
        from scipy.stats import chi2 as _chi2
        keep = np.where(mask)[0]
        X = honest[:, keep].astype(float)
        design = np.column_stack([np.ones(len(honest)), scores])
        bx, *_ = np.linalg.lstsq(design, X, rcond=None)
        Xr = X - design @ bx
        y = case.astype(float)
        by, *_ = np.linalg.lstsq(design, y, rcond=None)
        yr = y - design @ by
        xc = Xr - Xr.mean(0); yc = yr - yr.mean()
        den = np.sqrt((xc ** 2).sum(0) * (yc ** 2).sum()); den[den == 0] = np.inf
        r = np.clip((xc * yc[:, None]).sum(0) / den, -1 + 1e-12, 1 - 1e-12)
        pv = _chi2.sf(r ** 2 * len(honest), 1)
        hit = np.isin(keep, causal)
        rows.append({"lambda_gc": float(genomic_control_lambda(pv)),
                     "min_p_causal": float(np.min(pv[hit])) if hit.any() else None})
    return rows


def _boot_ci(diff):
    """10k-resample bootstrap 95% CI of a paired mean. The bootstrap RNG seed is
    fixed (7) so repeated runs of the same summary give identical intervals."""
    rng = np.random.default_rng(7)
    boot = diff[rng.integers(0, len(diff), size=(N_BOOT, len(diff)))].mean(1)
    return [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]


def summarize(res):
    """Paired gaps (vs all AF baselines, vs the conditioned single release,
    adaptive vs random) with bootstrap CIs."""
    for coh, byseed in res["cells"].items():
        seeds = list(byseed)
        af = np.array([byseed[s]["af_auc"] for s in seeds])
        afp = np.array([byseed[s]["af_pooled_auc"] for s in seeds])
        afc = np.array([byseed[s]["af_conditioned_auc"] for s in seeds])
        spc = np.array([byseed[s]["spec_conditioned_auc"] for s in seeds])
        summ = {"af_ref_mean": float(af.mean()), "af_ref_std": float(af.std()),
                "af_pooled_mean": float(afp.mean()), "af_pooled_std": float(afp.std()),
                "af_cond_mean": float(np.nanmean(afc)), "af_cond_std": float(np.nanstd(afc)),
                "spec_cond_mean": float(spc.mean()), "spec_cond_std": float(spc.std()),
                "strategies": {}}
        # the strongest passive baseline each strategy's gap is measured against
        vs_base = np.maximum(afc, spc)
        for st in STRATEGIES:
            summ["strategies"][st] = {}
            for q in QS:
                def pick(sd, key, q=q, st=st):
                    row = byseed[sd]["strategies"][st]
                    return row[q][key] if q in row else row[str(q)][key]
                auc = np.array([pick(s0, "membership_auc") for s0 in seeds])
                vis = np.array([pick(s0, "visible_share") for s0 in seeds])
                te = np.array([pick(s0, "target_err_outside") for s0 in seeds])
                ti = np.array([pick(s0, "target_err_inside") for s0 in seeds])
                rv = np.array([pick(s0, "rank_visible") for s0 in seeds])
                pm = np.array([pick(s0, "membership_perm_p") for s0 in seeds])
                se = np.array([pick(s0, "scatter_recon_err") for s0 in seeds])
                to = np.array([pick(s0, "target_resid_outside") for s0 in seeds])
                t01 = np.array([pick(s0, "membership_tpr_fpr01") for s0 in seeds])
                t05 = np.array([pick(s0, "membership_tpr_fpr05") for s0 in seeds])
                summ["strategies"][st][str(q)] = {
                    "membership_auc": float(auc.mean()), "membership_auc_std": float(auc.std()),
                    "membership_tpr_fpr01": float(t01.mean()),
                    "membership_tpr_fpr05": float(t05.mean()),
                    "gap_vs_af": float((auc - af).mean()),
                    "gap_vs_af_ci95": _boot_ci(auc - af),
                    "gap_vs_af_cond": float((auc - afc).mean()),
                    "gap_vs_af_cond_ci95": _boot_ci(auc - afc),
                    "gap_vs_strongest_passive": float((auc - vs_base).mean()),
                    "gap_vs_strongest_passive_ci95": _boot_ci(auc - vs_base),
                    "perm_p_max": float(pm.max()),
                    "visible_share": float(vis.mean()), "visible_share_std": float(vis.std()),
                    "scatter_recon_err": float(se.mean()), "scatter_recon_err_std": float(se.std()),
                    "target_err_outside": float(te.mean()),
                    "target_resid_outside": float(to.mean()),
                    "target_err_inside": float(ti.mean()),
                    "rank_visible": float(rv.mean()),
                }
        # strongest-passive comparators at 5% FPR (the single release and its
        # conditioned form both sit near spec_cond; report the passive release)
        summ["baseline_tpr_fpr05_single_release"] = float(
            np.mean([byseed[s]["strategies"]["passive"]["1"]["membership_tpr_fpr05"]
                     for s in seeds]))
        summ["spec_cond_tpr05_mean"] = float(np.mean(
            [byseed[s]["spec_conditioned_tpr05"] for s in seeds]))
        summ["no_release_controls"] = {
            tag: float(np.mean([byseed[s]["no_release_controls"][tag] for s in seeds]))
            for tag in byseed[seeds[0]]["no_release_controls"]}
        summ["adaptive_minus_random"] = {}
        for q in QS:
            a = np.array([byseed[s]["strategies"]["adaptive"][str(q)]["membership_auc"] for s in seeds])
            r = np.array([byseed[s]["strategies"]["random"][str(q)]["membership_auc"] for s in seeds])
            summ["adaptive_minus_random"][str(q)] = {
                "mean": float((a - r).mean()), "ci95": _boot_ci(a - r),
                "pos_splits": int((a > r).sum())}
        res.setdefault("summary", {})[coh] = summ
    return res


def run(chrom=22, seeds=5, ld_window=1500, max_snps=20000, caps=None):
    caps = caps or CAPS
    res = {"chrom": chrom, "seeds": seeds, "k_top": K_TOP, "block_len": BLOCK_LEN,
           "n_candidates": N_CANDIDATES, "n_perm": N_PERM, "qs": list(QS),
           "caps": caps, "cells": {}, "witness": {}, "power": {}, "summary": {}}
    out_name = f"s61_spectral_tomography_chr{chrom}.json"
    Ge, _, _ = s48.load_for(("EUR",), chrom, ld_window, max_snps)
    cohorts = {"within:EUR": [Ge, np.array(["EUR"] * len(Ge))]}
    G3, _, sup3, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    cohorts["super:3pop"] = [G3, sup3]
    G5, _, sup5, _ = s51._sup_label(SUPERPOP, chrom, ld_window, max_snps)
    cohorts["within:EAS"] = [G5[sup5 == "EAS"], sup5[sup5 == "EAS"]]
    for name, (G, groups) in cohorts.items():
        res.setdefault("marker_count", {})[name] = int(G.shape[1])
        for s in range(seeds):
            res["witness"][f"{name}:seed{s}"] = witness(G, groups, seed=s, caps=caps)
            res["cells"].setdefault(name, {})[f"seed{s}"] = cell(G, groups, seed=s, caps=caps)
            roles = split_roles(groups, caps, s)
            res["power"][f"{name}:seed{s}"] = gwaspower(G[roles["honest"]], s)
            _dumped(res, out_name)
    summarize(res)
    for name, summ in res["summary"].items():
        print("### %s  AF %.3f  AFpool %.3f  AFcond %.3f  specCond %.3f" % (
            name, summ["af_ref_mean"], summ["af_pooled_mean"], summ["af_cond_mean"],
            summ["spec_cond_mean"]), flush=True)
        for st in STRATEGIES:
            print("  %-8s " % st + " | ".join(
                "q%d auc%.3f gapStrong%+.3f vis%.3f" % (
                    q, summ["strategies"][st][str(q)]["membership_auc"],
                    summ["strategies"][st][str(q)]["gap_vs_strongest_passive"],
                    summ["strategies"][st][str(q)]["visible_share"]) for q in QS), flush=True)
    return res


def _dumped(res, name="s61_spectral_tomography.json"):
    """Write the JSON with NaN/Inf replaced by null so it stays strict JSON."""
    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [clean(v) for v in o]
        if isinstance(o, float) and not np.isfinite(o):
            return None
        return o
    with open(os.path.join(OUT, name), "w") as f:
        json.dump(clean(res), f, indent=2)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    payload = run(chrom=a.chrom)
    payload["elapsed_s"] = round(time.time() - t0, 1)
    _dumped(payload, f"s61_spectral_tomography_chr{a.chrom}.json")
    print("done", payload["elapsed_s"], "s", flush=True)
