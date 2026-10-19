# Round-12 Review (run_20261007_3a710cda)

**Reviewer**: pi (artifact-level verification). Read the four paper edits, the new/changed code
(s44, s46, s38 fnr_mode, enforce.py/s39 native injection, s37 nsp=16000), all six round-12
JSONs (incl. shard6 and the merged s37), the Round-12 results/tracker/plan sections.
**Date**: 2026-10-07

## The four adversarial questions — answered

**(2) Native-code test genuine, not mislabeled — YES.** `enforce.child_main` injects
`ctypes.CDLL(None)` under `INJECT_NATIVE=1` after the hook is installed; the hook raises
`SandboxViolation` on the `ctypes.dlopen`/`ctypes.dlsym` audit events, and the ONLY code path
that can set `native_egress_blocked=True` is that catch. JSON: 3/3 `native_egress_blocked=True`
with `reason='ctypes.dlopen'` (the hook's own event name — proof the hook fired, not an
unrelated failure) and `child_failed=False` (the child ran to completion; an unrelated crash
would have surfaced as `child_failed=True` and never produced a report). Paper's "denied on
ctypes.dlopen/ctypes.dlsym (3/3 seeds), so the native loading route does not escape the hook"
matches. The natural-egress/bypass matrix values are unchanged from round 11 ✓.

**(3) s46 uncoarsened comparator on the SAME object — YES.** `build_bundle2` places the
full-precision `rel["af_by_group"]` in the bundle; `features2(coarsen_site=False)` computes both
the single comparator and the joint-attacker features from exactly that object, with the same
`m3` (call-rate) features passed in (imported verbatim from `s40_bundle.call_rate_full`).
`single_coarse` (0.94707) reproduces s40's max_single (0.94707) — the "+0.026 against the
coarsened single object of the first manifest" claim is exact. Both comparisons are explicitly
labeled in the paper ✓.

**(4) s44 real/cohort genuine, controls at chance — YES.** Real arm: `real_multi.get_real_panel_chr`
(1000G chr22, 2241 SNPs), 8 splits via `three_way_split` (verified disjoint per super-population:
`idx[:per_p] / [per_p:2per_p] / [2per_p:3per_p]`). Cohort arm: `semi_synthetic.generate_cohort`
with the IDENTICAL call to s37's cached nsp2000_m8000 cell (same chrom/rho/M/seed), so the s37
LD/relatedness/eig controls at that cell cover it. Controls at chance: real 0.5021/0.5003,
cohort 0.4974/0.4990 ✓. The rule-3 floor is a legitimate confidence correction (FNR ≤ 3/n_pos,
the rule-of-three upper bound on a zero observed-FN count) — and the paper honestly reports that
the real-chr22 high-AUROC channels stay ordering-only (their ROCs hit TPR=1, so the strict check
fails regardless of floor; verified False at both δ).

**(1) Number match — all round-12 paper numbers verified against their JSON keys:**
native 3/3 ✓; s44 A* 0.27/0.35 real, 0.14/0.04 cohort, controls ≤0.07, "(8 and 1 splits)" ✓;
s46 0.973/0.964/+0.008/5-of-5/+0.026 ✓; s37 nsp=16000: 48,000 study genomes (n_site=48000=3×16000) ✓,
≈19× the 1000G arm (48000/2504=19.2) ✓, gap −0.0264 → "−0.03 to −0.10" range ✓, "narrows roughly
threefold before plateauing" ✓ (−0.096→−0.030 ≈ 3.2×, then 0.026–0.040), in-sample gap 0.0029 ✓,
controls pass (LD, relatedness 63.7 vs 62.5, eig10 0.213 vs 0.220) ✓; s38 resolution values
unchanged (paper's 5.0/4.1/1.1 still match) ✓; compile 16 pp ✓.

## Score: 8.5 / 10

## Remaining weaknesses (ranked)

**W12-1 (minor — paper, stale count).** §scale: "pairwise relatedness sits within ∼10% of the
donor library at **five of six** cells" — the delivered cohort grid now has **seven** cells
(nsp16000 added; its relatedness 63.7 vs 62.5 = +1.9%, within 10%). Should read "six of seven
cells". The rest of the paragraph (gap range, Limitations "up to 16000") was updated; this count
was missed.
**Minimum fix**: "five of six" → "six of seven".

**W12-2 (minor — docs).** `EXPERIMENT_RESULTS.md` B3: "LOO per-site AF drop .076" — that is
s40's figure. In s46 the LOO per-site-AF drop is 0.973 − 0.740 = **0.233** (manifest 2 has no
pooled AF, so dropping per-site AF leaves only the call-rate channel). The paper does not quote
the s46 LOO ✓, so the paper is clean.
**Minimum fix**: ".076" → ".233" (or label it as the s40 value).

**W12-3 (trivial — docs).** Results doc B4 is labeled "(S47)" but the tracker's round-12 table
has no S47 row (the native test was folded into s39/enforce.py; the table has S44/S45/S46).
**Minimum fix**: label B4 as S39/enforce.py, or add the S47 row.

**W12-4 (minor — paper, presentational).** The nsp=16000 headline cell runs at M=3000 — sparser
than the other cohort cells (6000–10000). Disclosed in the results doc/tracker/plan, but the
paper quotes "48,000 study genomes, ≈19×" for this cell without its marker count (consistent
with the paper's existing no-per-cell-M disclosure, but the headline cell deserves the clause).
**Minimum fix**: one clause — "at M=3000 markers" — after "nsp=16000".

**W12-5 (nit — paper, causal overstatement).** "Replacing the false-negative resolution with a
confidence-corrected floor (rule of three, 3/n_pos) **makes** the f-DP bound valid on the
cohort's every channel" — the resolution-floor bound was already valid on the cohort's every
channel (cohort `valid_res=True` everywhere; the floor never binds there — the cohort's ROCs
never reach TPR=1). The rule-3 floor's real role is a more conservative ε value; it does not
cause the cohort's validity. One-word fix: "admits a valid f-DP bound on the cohort's every
channel (with both floors)".

**Observation (unchanged, no action)**: EXPERIMENT_PLAN.md's removed B3 promises (documented in
the audit).

## Verified unchanged (regression)

Frontier table (0.590/0.964/0.588), k-ladder (0.53/0.65/0.74), all round-11 numbers (s38 eps
5.0/4.1/1.1 and controls ≤0.04/≤0.24; s39/s40/s41/s42/s37) — unchanged and matching. The s39
bypass-matrix values identical to round 11 (library leak 1.0, natural-egress blocked 1.0,
sanctioned blocked 1.0, chain 1.0, governed 0.624/0.509). Compile 16 pp, 5 overfull hboxes
(cosmetic).

## Verdict

**Almost — one one-word fix away.** All four round-12 blocks are genuine and verified: the
native-code boundary test is a real denial (not a mislabeled failure), the uncoarsened
comparator is computed on the same full-precision per-site AF as the bundle, the s44
real/cohort re-expression is real data with at-chance controls and an honest ordering-only
report on real chr22, and the nsp=16000 cell passes its controls at 48,000 study genomes. The
residuals are one stale count in the paper (W12-1) plus three doc-level nits. Apply W12-1
(one word) and W12-2..5 (one line each), recompile, and the next round is a pure regression
check.
