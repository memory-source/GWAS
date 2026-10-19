# Round-13 Review — regression check of W12-1..W12-5 (run_20261007_3a710cda)

**Reviewer**: pi (artifact-level verification; each fix grepped/verified in place, compile log
checked, key round-11/12 numbers re-grepped for regressions).
**Date**: 2026-10-07

## W12-1..W12-5 scorecard

| # | Fix | Status | Evidence |
|---|---|---|---|
| W12-1 | "five of six cells" → "six of seven cells" | **FIXED** | Paper: "pairwise relatedness sits within ∼10% of the donor library at **six of seven** cells (the smallest cell, nsp=150, is +11%)" — matches the 7-cell delivered grid (nsp16000 relatedness 63.7 vs 62.5 = +1.9%, within 10%). |
| W12-2 | Results-doc B3 LOO ".076" → ".233" | **FIXED** | Now "LOO per-site AF drop **.233**" — matches s46 (bundle2 0.9726 − loo af_site 0.7395 = 0.233). The ".076" remaining in the round-11 B3 section is s40's correct value in its own context. |
| W12-3 | B4 header "(S47)" numbering | **FIXED** | Now "## B4 — Native-code egress boundary test (S47 → merged into S39/enforce.py)" — the numbering is explained. |
| W12-4 | Paper "at M=3000 markers" after nsp=16000 | **FIXED** | "the largest delivered cell is nsp=16000 at M=3000 markers (48,000 study genomes, ≈19× the 1000G arm…)" — matches the JSON (nsp16000_m3000, n_site=48000, M=3000). |
| W12-5 | "makes the f-DP bound valid" causal overstatement | **FIXED** | Now "A confidence-corrected false-negative floor (rule of three, 3/n_pos) **admits** a valid f-DP bound on the cohort's every channel" — accurate (cohort valid_r3=True on all channels; the resolution-floor bound was already valid there, so "admits" is the right verb). Sentence flows correctly into the ordering-only real-chr22 caveat. |

## Regression check

Compile: 16 pp, **0 errors** (5 overfull hboxes, cosmetic). The round-12 numbers are all still
in place (s44 0.2714/0.350, cohort 0.14/0.04, s46 0.973/0.964/+0.008, 48,000 study genomes,
≈19×, ctypes.CDLL(None) denied 3/3) and the round-11 numbers unchanged (s38 ε 5.0/4.1/1.1,
0.55/4.1, A* direction, s39 0.624/0.509). No new edits beyond the five fixes.

## Score: 8.5 / 10

## Remaining weaknesses

None found. The remaining items are the paper's own disclosed limitations, each already carrying
its caveat in the text: all membership-level arms simulated or semi-synthetic (real 1000G with
synthetic splits; the biobank cohort is mosaic-copying); the missingness channel
constructed-by-design with the ρ=0 null exactly at chance (0.502) as the control; 3 seeds/cell
in s37/s41/s42 (tight stds on the headline deltas); ε_lb presented as ordering information
wherever a zero-observed-FN corner exists (flags disclosed). Minimum fixes: none required.

## Verdict

**YES — ready for submission** as a measurement paper. Every defect found across the four
review rounds (W1–W8 in round 11, R1–R5 in loop 2, W12-1–W12-5 in round 12) is now fixed and
verified in the artifacts; every round-11/12 number matches its JSON key; the native-code
boundary test, the uncoarsened comparator, the real/semi-synthetic per-individual re-expression,
and the biobank-corner cell are all genuine. Recommend: no further review rounds needed; any
future change should be re-verified against the JSONs as in this loop.
