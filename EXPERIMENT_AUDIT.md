# Experiment Audit Report

**Date**: 2026-10-14
**Auditor**: pi-review / glm-flash, ultra reasoning (cross-model, read-only)
**Project**: Audit-Governed Multi-Site GWAS (reveal-leakage measurement + governance)
**Round**: 17 — s59 GPU-AUROC tie disclosure, s60 monotone share + per-seed crossing,
s56 fixed-reference rewrite, s54/s55/s37/s42 text corrections
**Trace**: thread `2601f76071434a2f8e8403eb2dfe1126`

## Overall Verdict: WARN

## Integrity Status: warn

No fraud found in the measurement pipeline: ground truth is genuine (real 1000G
3-way disjoint splits), splits/orientation are clean, no output-derived
normalization, no dead code, and headline numbers reproduce to the printed digit.
The WARN is for paper/doc statements that contradicted their own JSONs and for a
tautological result dressed as a check — all corrected in this round.

## Checks

### A. Ground Truth Provenance: PASS
Membership is a real 3-way disjoint per-group split (`s48_real_structure.py:44-58`,
`s56_real_poolmax.py:44-60`, `s60_share_crossing.py:31-53`). s54 reproduces s48's
AUROCs to all printed digits. s37/s42 proxy arms are explicitly labeled.

### B. Score Normalization: WARN
s60's denominator (total squared norm about the study mean) is an input statistic,
not model output. But the cumulative squared-projection share is monotone **by
construction** (cumsum of non-negative increments); the paper's original sentence
cited that monotonicity as if it certified a gap. Reworded: it only rules out a
share that rises by construction.

### C. Result File Existence: WARN
Every headline number pairs to a live JSON key, but two quoted statements were
false against their JSONs and were fixed:
- s42 (`paper/main.tex`): "differ by <=0.01 / both bounded above by candonly" —
  actual diff 0.0212 at rho=0.25 and table 0.893 > candonly 0.888 at rho=0.75. Now
  "track within ~0.02, marginally above at rho=0.75".
- s56: "up to sd 0.26" — actual max 0.41 (chr1 EAS:n350). Corrected in paper +
  results doc.
- "rather than all eighteen" — the s48 ladder has 19 cohorts. Corrected.

### D. Dead Code Detection: PASS
s59/s60 functions are all called and reach JSON. One harmless dead assignment
(s60 `crossing.mean` overwritten) leaves the JSON intact.

### E. Scope Assessment: WARN
- s60 crossing: per-seed instability now disclosed (seed 0 triple-crosses, seed 4
  does not cross; mean ~341 over 3 single-crossing seeds).
- s56: reference held at 150 (capped 147-148 only at the pool boundary); EUR/EAS
  n350 cells have only 6 non-members, now stated in the paper; the AFR n350 tuple
  quoted is coherent.
- s54: the permutation null covers the five tested cohorts, not the full ladder.

### F. Evaluation Type
real_gt for s54/s55/s56/s60; simulation_only for s59; synthetic_proxy (labeled)
for s37/s42. No mislabeling.

## Action Items (all applied this round)
- Fix the s42 "<=0.01 / bounded above" claim (paper + tracker note).
- Fix the s56 "sd up to 0.26" claim (paper + results doc) and state the n_non=6
  degeneracy for the EUR/EAS n350 cells.
- Reword the s60 monotone-share sentence to a non-tautological statement and
  disclose the per-seed crossing instability in the paper.
- Add the s59 tie disclosure to the paper Limitations (was tracker-only) and
  correct the "eighteen" count; fix the stale round-15b s56 tracker row.

## Claim Impact
- Structure-zero PC>AF bound (AFR n350): supported, scoped.
- s56 fixed-reference design: supported.
- s60 monotone share: needs_qualifier (now stated as construction, not evidence).
- s59 GPU tie: supported (bounded disclosure; archived cells tie-free).
- s42 candidate-only match: supported after the correction.
