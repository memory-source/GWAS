Round 2/4 — re-review after implementing all four round-1 weaknesses.

Read the changed artifacts yourself; do not take my word for what changed or whether it
worked:
- `code/s72_printed_matched.py` + `results/raw/s72_printed_matched_chr{22,1}.json` (NEW, W1)
- `code/s67_joint_passive.py` + `results/raw/s67_joint_passive_chr{22,1}.json` (W4 + TPR bar)
- `code/s70_known_fraction.py` + `results/raw/s70_known_fraction_chr{22,1}.json` (W4)
- `code/s71_protection_cost.py` + `results/raw/s71_protection_cost_chr{22,1}.json` (W2)
- `code/s69_reliable_dirs.py` (W3)
- `refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md`, `CLAIMS_FROM_RESULTS_ARB_20261017.md`
- Reviewer memory: `review-stage/REVIEWER_MEMORY.md`

W1 (matched printed comparison): s72 runs q = 8 joint re-fits per cohort x split through
real PLINK 2.0 `--pca 10 allele-wts approx --seed 1000+t`, reads the printed weights and
eigenvalues back, forms `Y_t = V_t Λ_t − K_t V_t` in PLINK's scatter/M normalization
with the partner's `K_t` known exactly, and compares the q8 stack against the four
passive constructions formed from the SAME printed outputs. The passive adder does not
hold `K_t`, so it forms the conditioned single release and fuses over q. Result: q8
beats the printed passive bar by +0.0070 EUR chr22 / +0.0126 EUR chr1 (positive CI) and
+0.0096 chr22 / +0.0101 chr1 3pop (CI spanning zero); `oracle_q ≈ oracle_single`, so
handing the passive adversary all q outputs does not close the gap. The matched
cross-population TPR gap is negative, so J9 is labeled an AUC result.

W2: `fp_rate` renamed `off_causal_fraction`; the "residual structure" attribution
removed; "usable" replaced with a low-power stability statement.

W3: J5 reworded to "heuristic retained rank under simulated quantization (offline float64
reference)"; the results doc now states the quantization is decimal rounding of unit-norm
eigenvector entries and is NOT PLINK's six-significant-digit print.

W4: s67 and s70 persist `per_split` rows; both JSONs carry a `ci_definition` string
labeling the intervals as conditional split uncertainty over one fixed panel.

Also incorporated your round-1 comparator note: the TPR reading now uses a
metric-specific TPR oracle over the four candidates, reported alongside the AUC oracle
(`tpr_oracle_bar`, `gap_vs_tpr_oracle`).

Please re-score and re-assess. Are the four weaknesses addressed? Same format: Score,
Verdict, Remaining Weaknesses, Minimum Fixes. Include a Memory Update at the end.
