# Round-15 Review, Loop 2 (re-review after the R15-1 6/10 verdict) — run_20261010_a1f2c3d4

**Reviewer**: pi (artifact-level verification; adversarial re-read of code + JSON + paper).
**Date**: 2026-10-10.
**Inputs read in full**: `paper/main.tex` (21 pp, fresh compile 0 errors), `code/s51_real_lowstruct.py`,
`s52_haplo_second.py`, `s53_eigbasis.py`, `s38_privacy_metric.py`, `s37_fitsens.py`, plus the shared
helpers they depend on (`s48_real_structure.cohort_cell/_orient_auroc`, `s32_phased_ibd.plant_phased_ibd`),
the eight result JSONs, `CLAIMS_FROM_RESULTS.md`, `refine-logs/EXPERIMENT_RESULTS.md`, tracker rows,
and `s10_paired.json` for the minor-fix check.

## 1. The four findings from R15-1 — verified fixed (in the paper)

1. **eps validation circular → FIXED.** `validate_tradeoff` now applies the *same* FNR floor as
   `eps_lower_bound` (s38:150–171); the paper says explicitly "the `valid` flag checks consistency
   with that convention, not statistical validity" (main.tex:706), names the default 1/n_pos floor
   and the rule-of-three 3/n_pos reading, and states the quoted ε values use the default floor.
   The inversion example is real and stated: af_full AUC .8446 / ε 4.120 vs af_2dp .8396 / ε 4.190
   (JSON `eps_order_inversion`). The confidence-upgrade wording is gone. The `valid` flags are no
   longer used to distinguish real vs cohort (see 4).
2. **"Pure total-size effect" → FIXED in the paper.** Heading is now "Fixed super-population
   proportions, growing size"; the text says "jointly increasing the study, reference and
   evaluation panels produces a decreasing split-mean gap … this is a joint size effect rather
   than an isolated study-size effect", and the crossing locations (n_study≈260–400) are labeled
   "interpolations of the five-split means". Both interpolated crossings recompute exactly
   (chr1 ≈258 between per80/per120; chr22 ≈385 between per120/per160).
3. **Pooled arm naming → FIXED.** The paper explains the definitional zero ("pools each
   super-population into one group so measured between-group structure is zero by construction"),
   qualifies the cap ("n_study=167–220 … a property of that split design, not an absolute 1000G
   ceiling" — matches n_pool 503/504/661), and keeps the honest "zero-structure real AF-first
   cohort remains unobserved inside 1000G".
4. **S44 staleness → FIXED.** `s44_individual_real.json` re-run (Oct 10 00:58, after the s38 helper
   change 23:11); every `*_valid` flag in s38/s44 is True; the paper now reads ε uniformly
   ("descriptive and floored under a stated convention … on every channel"). The surviving
   Limitations clause "(ordering-only on the real high-AUROC channels)" is a resolution statement,
   not the obsolete validity distinction — correct as written.
5. **Minors → FIXED.** S10 "0.941 ± 0.012" (JSON sd .012465 ✓); msprime swapped null 0.518 ✓;
   "the full 48000-row fit was not run" (no GB claim) ✓; s53 docstring renamed "deterministic
   cosine directions (… not exactly an orthonormal DCT)" and the paper's "0.48–0.52" range now
   covers DCT-EAS 0.484 ✓; s37 docstring corrected ✓.

## 2. The five limitations — implemented as claimed, numbers trace 1:1

- **L1 (s51)**: pooled within-superpopulation arm (chr22 +0.075/+0.051/+0.044; chr1
  +0.110/+0.165/+0.048; controls 0.465–0.526) and the fixed-composition ladder (chr22
  +0.120/+0.042/+0.007/−0.026 at Fst 0.0873–0.0911; chr1 +0.058/+0.002/−0.010/−0.018 at
  Fst≈0.103–0.104). Every quoted number matches its JSON key. The s48 `cohort_cell` machinery uses
  a held-out-half orientation with a shared calibration mask; ref/eigonly controls at chance.
- **L2 (s38)**: floor convention + inversion as above; A*/ε table matches (PC .610/5.00; AF full
  .546/4.12; AF 1dp .187/1.06; controls ≤.067/≤.244).
- **L3 (s52)**: two additional leave-self-out statistics; real arm at chance (chr22 .488/.496,
  chr1 .497/.497); planted control now *recovers* (chr22 .798/.797, chr1 .776/.788) with si-only
  planting and seg=1500 ≈ 31% of the panel — the planted AUROC ≈ .80 is the expected ceiling for a
  60%-planted member arm (0.6·1.0 + 0.4·0.5), i.e. the control detects planted members essentially
  perfectly. Self-block removal is index-correct in the planted arm (donor = member 0, unmodified,
  remains in the release).
- **L4 (s53)**: identity / deterministic-cosine / λ-weighted fixed bases all at chance
  (.484–.517) while study-fitted loadings give .81–.89; the λ-weighted score divides by the
  released eigenvalues, which *is* the tested control (documented).
- **L5 (s37_fitsens)**: gap −0.030 (sd .008) at every n_fit incl. the full 24000-row fit;
  −0.027/−0.027/−0.026 at nsp16000; "was not run" stated honestly; af_only constant across n_fit
  as expected.

## 3. Score: 7.5 / 10 — ALMOST

The substance holds up under adversarial reading: every new quoted number traces to a JSON key,
the fixed prior findings stay fixed, and the two controls that had failed earlier (s52 planted,
s44 validity) now genuinely work. What remains is propagation/hygiene, not science: the s53
robustness result is reported only in the abstract, and the internal documents still carry the
exact phrasings the previous review rejected.

## 4. Remaining weaknesses (ranked; minimum fixes)

- **W1 (minor–moderate; the one real paper gap).** The s53 fixed-basis result exists **only in the
  abstract** (main.tex:73–75). The mechanism paragraph ("The leak needs the loadings…",
  :410–416), the Discussion's governed-object implication (:1169–1172), and the §governance
  per-individual paragraph still describe the eigenvalue-only negative as the *random-basis*
  ablation without the strengthening evidence. Nothing is false — but the paper's best support for
  its punchline (withhold loadings; eigenvalue-only is inert) lives where skimmers don't read.
  **Fix**: one clause in the mechanism paragraph and one in the Discussion, e.g. "fixed
  data-independent bases, including an eigenvalue-weighted one reusing the released spectrum, are
  also at chance (s53_eigbasis.py, 0.48–0.52)". Two sentences, no new experiments.
- **W2 (minor; stale internal docs).** `EXPERIMENT_RESULTS.md` L1 bullet still says "the real
  PC-vs-AF switch is a **pure total-size effect** at fixed composition, crossing at n≈400" — the
  exact phrasing R15-1 rejected (paper now: joint size effect, crossings ≈260–400); its L5 bullet
  retains "(full Gram >24 GB)" which the paper corrected to "was not run". Tracker L4 row says
  ".49–.52" but DCT-EAS is 0.484. **Fix**: three one-line doc edits.
- **W3 (minor; claims doc lag).** `CLAIMS_FROM_RESULTS.md` (last updated round 25) does not mention
  s51/s52/s53 anywhere. Its "do not support" entries remain accurate in substance, but C6 is now
  directly supported by the fixed-composition ladder and the evidence list omits three JSONs.
  **Fix**: add s51/s52/s53 to the Evidence line + one sentence under C6 and under the
  "no real low-structure large-n cell" entry (pooled arm shows PC>AF at definitional-zero
  structure up to n=220; zero-structure AF-first at size remains unobserved).
- **W4 (nit; wording).** "Two further, **independent** leave-self-out haplotype statistics" — they
  are two functional forms (max vs count>0.78) of the *same* per-haplotype IBS-share matrix, hence
  strongly dependent. The claim that matters (the negative is not specific to the LSR estimator)
  holds. **Fix**: "independent" → "additional".
- **W5 (nit; stale code comment).** s38 `run()` comment "(2dp has the higher AUROC but the smaller
  eps at every delta)" is doubly false (af_full has the higher AUROC *and* the lower ε); the JSON
  note beside it is correct. **Fix**: one-line comment edit (this loop fixed two identical-class
  comments in s10 at C12).
- **W6 (cosmetic).** Two overfull hboxes of 76–89pt (the s52 paragraph and Limitations) are visible
  in the PDF; 21 pp, 0 errors, 0 undefined refs.

Optional (not required): one clause disclosing THR=0.78 in s52 as fixed a priori — immaterial in
practice because the threshold-free max_ibs statistic agrees.

## 5. What plainly holds (verified, not trusted)

The disjoint-negative design; the held-out-half orientation machinery in s48/s49/s51/s53; the
planted-control arithmetic (0.8 ceiling, don't flag near-1.0 or near-0.8 as suspicious); the
floor-convention consistency between `eps_lower_bound` and `validate_tradeoff`; the s51 ladder's
honest disclosure chain (joint effect, interpolation, split-design cap, zero-structure unobserved);
the s37 fit-fraction invariance including a genuine full fit at 24000 rows. The manuscript states
exactly what its artifacts show at every position I checked.

**Verdict**: ALMOST. One consistency pass (W1's two sentences, W2/W3 doc edits, W4/W5 one-liners)
away from READY; no new experiments required.
