# Experiment Audit Report

**Date**: 2026-10-07
**Auditor**: pi-review / glm-flash, ultra reasoning (cross-model, read-only)
**Project**: Audit-governed multi-site GWAS — round 12 (s44/s45/s46/s47 additions)

## Overall Verdict: PASS

## Integrity Status: warn

## Checks

### A. Ground Truth Provenance: PASS
`real_multi.three_way_split` (code/real_multi.py:118-129) gives non-overlapping
per-super-population slices; the s44 real arm (`code/s44_individual_real.py:83-91`)
uses `get_real_panel_chr(22)` on the real 1000G Phase-3 VCFs and the cohort arm
uses `semi_synthetic.generate_cohort` (real phased donors) with labels assigned
after generation. `s37_realsynth_scale.py` nsp=16000 gives `n_site=48000` ✓.
`s40_bundle.call_rate_full` is a constructed released-object feature, disclosed,
never a label source. No label derives from model output.

### B. Score Normalization: PASS
No metric divides by model-output statistics. The `rule3` floor (`3/n_pos`) caps
the second branch of the f-DP bound, replacing round-11's unsound `1e-12` clip
(which produced ε≈24); it lowers ε_lb — the conservative direction for a
lower-bound claim. Both floors and `*_valid` flags are reported side by side with
AUROC/TPR@1%/KS/ρ/LR_max.

### C. Result File Existence / Number Match: WARN (one stale tracker row)
All five round-12 artifacts exist, are keyed as cited, and were independently