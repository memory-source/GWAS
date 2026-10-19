# Round-22 auto-review-loop — Round 2 inputs (post-fix)

Round 1 verdict was 6/10 "No" with four weaknesses. All four were implemented.

## W1 (J9) — the matched printed-output comparison (most important)

New script `code/s72_printed_matched.py`; raw `results/raw/s72_printed_matched_chr{22,1}.json`.

Per cohort x split (5 seeds), q = 8 joint re-fits through real PLINK 2.0
(`--pca 10 allele-wts approx --seed 1000+t`), read the printed allele weights and
eigenvalues back, form `Y_t = V_t Λ_t − K_t V_t` from the READ-BACK V/Λ in PLINK's
scatter/M normalization with the partner's own scatter `K_t` known exactly, and compare
the q8 stack against the passive constructions built from the SAME printed outputs.
The passive adversary does not hold `K_t`, so per release it forms the conditioned single
release `Z_k V_t` removal and fuses over q. Two matched bars: `oracle_single` (best of
the four constructions of the t=0 output) and `oracle_q` (best of the four fused over all
q printed outputs).

Results (AUC gap of q8 over each bar, 10k bootstrap 95% CI, positive split count):

| cohort | chr | oracle_single | oracle_q | q8 | gap vs oracle_single | gap vs oracle_q | pos |
|---|---:|---:|---:|---:|---:|---:|---:|
| within:EUR | 22 | 0.992 | 0.992 | 0.999 | +0.0070 [+0.0006,+0.0133] | +0.0070 [+0.0006,+0.0133] | 4/5 |
| super:3pop | 22 | 0.840 | 0.840 | 0.850 | +0.0096 [-0.0072,+0.0298] | +0.0097 [-0.0071,+0.0298] | 4/5 |
| within:EUR | 1 | 0.984 | 0.984 | 0.997 | +0.0126 [+0.0057,+0.0186] | +0.0125 [+0.0057,+0.0183] | 4/5 |
| super:3pop | 1 | 0.846 | 0.846 | 0.856 | +0.0101 [-0.0092,+0.0269] | +0.0101 [-0.0092,+0.0269] | 3/5 |

The passive bar does NOT strengthen when handed all q printed outputs (`oracle_q ≈
oracle_single`), so the increment persists in printed output. The matched cross-pop
TPR@5%FPR gap is NEGATIVE at this size (-0.014 chr22, -0.059 chr1); J9 is an AUC result.

## W2 (J7) — off-causal discovery fraction, not a clean false-positive rate

`code/s71_protection_cost.py` renamed the field `fp_rate` -> `off_causal_fraction` and
its counter `n_fp` -> `n_off_causal_discoveries`. The section docstring now states that
under LD some off-causal markers are proxies of causal ones, so a discovery there is not
a false positive. Results re-run: lambda/power/recovered-count/off-causal/agreement all
unchanged between the fixed and input-varying pass. The word "usable" was replaced by a
low-power stability statement at one study size. Raw:
`results/raw/s71_protection_cost_chr{22,1}.json`.

## W3 (J5) — the reliability floor is a heuristic count

`code/s69_reliable_dirs.py` docstring now states the perturbation `dY` is measured
against the float64 operand, which the attacker never sees, so the floor is an offline
construct for the measurement and the rank numbers are descriptive counts under one
plausible floor, not a certified usable-rank theorem. Results doc J5 heading reworded to
"the reliable-direction set as a function of export precision (heuristic)".

## W4 — per-seed measurements persisted

`code/s67_joint_passive.py` and `code/s70_known_fraction.py` now write `per_split` rows
in each cohort block. Re-ran chr22 and chr1. Raw:
`results/raw/s67_joint_passive_chr{22,1}.json`, `s70_known_fraction_chr{22,1}.json`.

## Updated narrative

`refine-logs/EXPERIMENT_RESULTS_ARB_20261017.md` (new J9 section, J4/J5/J7 rewording,
J9 added to "What these results support", printed-output TPR caveat added to "do not
support"); `CLAIMS_FROM_RESULTS_ARB_20261017.md` (center claim updated: the increment
now holds in printed output per the matched cell; J9 row flipped to yes for AUC, not
established for TPR; J5/J7 reworded); `refine-logs/EXPERIMENT_TRACKER.md` (S72 row,
D22b row).
