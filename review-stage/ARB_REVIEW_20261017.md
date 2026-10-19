# Round-22 auto-review-loop log (selection-input spectral tomography — fair comparison + interface)

Reviewer: codex / gpt-6-astra (read-only, cross-family vs deepseek-flash executor),
thread `01a149b1-6c3b-7f30-b275-e5c2e8c9230e`. MAX_ROUNDS = 4.
Artifacts: s67/s68/s69/s70/s71/s72 + `results/raw/*_chr{22,1}.json`,
`refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md`, `CLAIMS_FROM_RESULTS_ARB_20261017.md`.

## Round 1 — 6/10, "No"

Round-1 weaknesses (ranked):
1. The practical attack result lacked the decisive matched comparison — s67 measured the
   four-construction comparison in float64; no cell measured the accumulated attack
   against the passive comparator in real printed PLINK output.
2. J7's "false-positive" interpretation was unsupported and "usable" overstated; off-causal
   discoveries are LD-correlated with causal markers, not clean nulls.
3. J5 was stronger than its simulated-quantization experiment warranted.
4. The "raw results" files discarded the per-seed measurements, so the intervals could not
   be audited from the JSONs.

## Round 2 — 6/10, "No" (after s72, per-split, off-causal rename, s69 scope)

Genuine improvements confirmed (J5 scoped, per-split persistence, metric-specific TPR
oracle). One critical + one major + one medium + one lower remained:
1. **Critical**: s72 fitted the `fused_cal` logistic weights on CONDITIONED calibration
   features (already-member-subtracted), weakening the comparator. The reviewer restored
   s67's recipe (weights on UNCONDITIONED features) and found the corrected chr22 EUR gap
   +0.005429 [−0.000277,+0.011013] — the advertised positive interval did not survive.
2. **Major**: the all-q interpretation did not follow from the executed information model;
   `oracle_q` was a late-fusion diagnostic, not general processing.
3. **Medium**: s72 export paths omitted the cohort, so 3pop overwrote EUR.
4. **Lower**: the TPR split count was misstated (EAS 4/5 chr22, 3/5 chr1); "nominal level"
   and "individual QC is preserved" overstated.

## Round 3 — 7/10, "Almost" (after comparator repair, path tagging, wording)

The reviewer independently replayed chr1 EUR seed0 from eight printed exports, reproduced
AUROC 0.9968443 and conditioned-single 0.9828191, verified all 16 gap summaries and 160
exports, and confirmed no numerical defect. Remaining: text-only — headline cross-population
qualification, three stale TPR/count sentences, J7 nominal-error wording, tracker cohort
count.

## Round 4 — 7/10, "READY"

Final verdict: **READY for the bounded mechanism-and-measurement submission**. No remaining
weakness load-bearing for the narrowed claim. One nonblocking wording note (distinguish
split counts from intervals for mean effects), applied.

## Final claim state

- Float64: the accumulated q8 attack beats the metric-specific passive bar over four
  constructions in all three cohorts on both chromosomes, positive AUC CI (5/5 splits,
  chr22 EAS 4/5). Cross-population carries the large increment.
- Printed output (s72): positive-mean AUC gap in every cell (+0.0054 EUR chr22, +0.0082
  3pop chr22, +0.0100 EUR chr1, +0.0060 3pop chr1); conditional interval excludes zero only
  for chr1 EUR. The late-fusion passive bar does not strengthen with the q outputs. The
  matched TPR@5%FPR reading is not established.
- J5 is a heuristic retained-rank count under simulated quantization with an offline
  float64 reference.
- J7 is a summary-level stability statement at one study size; the off-causal discovery
  fraction is not a null-error estimate.
- J8 (independent subjects) not run.
