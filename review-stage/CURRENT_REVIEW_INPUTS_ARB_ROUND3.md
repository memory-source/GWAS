# Round-22 auto-review-loop — Round 3 inputs (post-Round-2 fixes)

Round 2 verdict was 6/10 "No" with four weaknesses. All four were implemented. W1's
corrected numbers reproduce the reviewer's independent computation exactly.

## W1 (critical) — s72 restored the s67 calibrated-comparator recipe

`s72_printed_matched.oracle_bar` now fits the `fused_cal` logistic weights on the
UNCONDITIONED calibration features (`_z(s_pooled)`, `_z(spec_uncond)`) and applies them to
the CONDITIONED channels, exactly as `s67_joint_passive.py` does. Re-ran chr22 + chr1.
The corrected chr22 EUR gap is **+0.0054, CI [−0.0003, +0.0110], 4/5** — matching the
reviewer's reconstructed +0.005429 [−0.000277, +0.011013]. Corrected chr22 3pop is
+0.0082 [−0.0084, +0.0287], also matching. Raw:
`results/raw/s72_printed_matched_chr{22,1}.json`.

Full corrected table (AUC gap of q8 over each bar):

| cohort | chr | oracle_single | oracle_q | q8 | gap vs single | gap vs q | pos |
|---|---:|---:|---:|---:|---:|---:|---:|
| within:EUR | 22 | 0.993 | 0.993 | 0.999 | +0.0054 [−0.0003,+0.0110] | +0.0055 [−0.0001,+0.0110] | 4/5 |
| super:3pop | 22 | 0.842 | 0.842 | 0.850 | +0.0082 [−0.0084,+0.0287] | +0.0082 [−0.0084,+0.0287] | 4/5 |
| within:EUR | 1 | 0.987 | 0.987 | 0.997 | +0.0100 [+0.0038,+0.0165] | +0.0100 [+0.0037,+0.0162] | 4/5 |
| super:3pop | 1 | 0.850 | 0.850 | 0.856 | +0.0060 [−0.0149,+0.0237] | +0.0060 [−0.0149,+0.0236] | 3/5 |

The claims file now states: J9 is PARTIAL — a positive-mean AUC gap in every cell,
but the conditional interval excludes zero only for chr1 EUR. The "positive CI on EUR
chr22" statement was removed.

## W2 (major) — all-q interpretation corrected

The s72 comment block now states the passive comparator is built under the SAME
information model as s67's `spec_cond` (it sees the released V/Λ and the operand Y_t,
which subtracts the known calibration members' projection), and `oracle_q` is labeled a
**late-fusion diagnostic over the four specified constructions**, not a general statement
about arbitrary processing of the transcript. The "outsider does not hold K_t" narrative
was removed.

## W3 (medium) — cohort-tagged export paths

`s72` export paths now include the cohort tag (`s72_c{chrom}_{cohort}_{seed}_{t}`), so the
EUR and 3pop exports no longer collide.

## W4 (lower) — counts and QC wording

- J1 TPR count corrected: "5/5 within-population TPR splits except chr22 EAS 4/5 and
  chr1 EAS 3/5".
- "sit at the nominal level" removed.
- QC statement now says precisely: the marginal per-SNP quantities (histogram, AF,
  variance, scatter) are exactly preserved, and individual-level statistics change
  slightly — max row-heterozygosity change 4.5e-3, max per-SNP variance change 2.2e-16 to
  6e-16 (machine precision).

Changed: `code/s72_printed_matched.py`, `refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md`,
`CLAIMS_FROM_RESULTS_ARB_20261017.md`, `refine-logs/EXPERIMENT_TRACKER.md`. Raw:
`results/raw/s72_printed_matched_chr{22,1}.json`.
