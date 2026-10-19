# Reviewer Memory — Selection-Input Spectral Tomography round (deep-search, 2026-10-16)

Backend: codex (gpt-6-astra, read-only, cross-family). This is a NEW research
direction round built on `research-plan/20261016_deep_search/DEEP_DIRECTIONS_ZH.md`,
separate from the passive-paper rounds.

## Threads
- codex deep-review loop: `01a14176-aaa4-7bf0-bf4c-9065eda1a1fe` (10 rounds, 4 → 7/10, verdict ready)
- pi-review loop: `54e81425ae86484ba6c4e751d772fe15` (2 rounds, 6.5 → 7/10, verdict almost)
- codex round-2 integrity re-audit: `01a14455-c02b-7c93-b103-24ed8aace24a` (FAIL → WARN → PASS)

## Artifacts
- Plan: `refine-logs/EXPERIMENT_PLAN_DEEP_20261016.md`
- Code: `code/s61_spectral_tomography.py`, `code/s62_sample_size.py`,
  `code/s63_plink_export.py`, `code/s64_saturation.py`,
  `code/s65_identity_ablation.py`, `code/s66_control_provenance.py`
- Raw: `results/raw/s61_spectral_tomography_chr22.json` / `_chr1.json`,
  `s62_sample_size_chr22.json`, `s63_plink_export_chr22.json`,
  `s64_saturation_chr22.json`, `s65_identity_ablation_chr22.json`,
  `s66_control_provenance_chr22.json`
- Results: `refine-logs/EXPERIMENT_RESULTS_DEEP_20261016.md`
- Audit: `EXPERIMENT_AUDIT_DEEP.md` / `.json` (PASS)
- Claims: `CLAIMS_FROM_RESULTS_DEEP_20261016.md` (verdict yes)
- Review log: `review-stage/DEEP_REVIEW_20261016.md`

## Substantive lessons (rounds 1–3)
- A conditioned single-release baseline must fix its orientation on the UNCONDITIONED
  release; fixing it after removing calibration members flips the sign per split.
- AUC on a ceiling-censored channel (attack at 1.000) understates the practical
  effect; report TPR@5%FPR alongside.
- A random-span control can be an artifact when the panel is MAF-filtered on a larger
  sample set than the cohort: globally-common/EAS-rare markers have cohort `mu == 0`
  and a scale clip gives out-of-`mu` carriers huge |z|; bound |z| for control rows.
- A control's conclusion must be shown under the control's own preprocessing choice,
  not only the pipeline's — span-fitting and scoring must share one transform.
- An ablation that supports a claim belongs in an artifact (JSON), not reviewer
  scratch.
