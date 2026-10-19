# Round-11 Review — Loop 3 (run_20261006_f68b0cd6)

**Reviewer**: pi (artifact-level verification; cumulative diff of `paper/main.tex` vs the
round-1 snapshot = 7 hunks, every one mapping to a reviewed fix; JSON/docs/tracker re-read).
**Date**: 2026-10-07

## R1–R5 scorecard (round-2 residuals)

| # | Residual | Status | Evidence |
|---|---|---|---|
| R1 | Paper quoted af_2dp's values (0.53/4.2) under the "full-precision AF" label | **FIXED** | Paper now: "versus $0.55$/$4.1$ for full-precision AF and $0.19$/$1.1$ at 1dp" — af_full 0.5456/4.120 ✓, af_1dp 0.1871/1.057 ✓. Remaining "0.53" hits are unrelated (k-ladder, eig-only range). |
| R2 | B5 CoV line garbled ("donor CoV 186/141") | **FIXED** | Now: "at nsp150 cohort CoV 186/141=1.32 vs donor 48/128=0.38, other cells' cohort CoV ~0.3". Verified: cohort 186.03/141.25=1.317 ✓, donor 48.12/127.74=0.377 ✓, other cells 0.281–0.362 → "~0.3" ✓ fair. |
| R3 | Tracker S38 stale ">0" shorthand | **FIXED** | Now "A*−(2AUC−1)<0 leaking / >0 controls (non-binormal)" — matches JSON. |
| R4 | s41 "at full marker density" (W6 remainder) | **FIXED** | Now "at uniform genome-wide density", parenthetical "uniformly strided from the full-density parse" retained — sentence internally consistent. |
| R5 | `eps_lb_*_valid=False` flags undisclosed | **FIXED** | B1 clause added: "For the four high-AUROC channels the `eps_lb_*_valid` flag is False because a single zero-observed-FN point makes the strict empirical trade-off infeasible for any finite ε; the 1/n_pos floor is the convention that resolves it, so those values are read as ordering information." Verified: exactly af_full/af_2dp/af_sitewise/pc_loadings carry valid=False (af_1dp + both controls True) — clause accurate. |

## Cumulative status

All eight round-1 weaknesses (W1–W8) and all five round-2 residuals (R1–R5) are now addressed
and verified in the artifacts. Regression: the paper diff vs round 1 contains only the reviewed
fix hunks; s39/s40/s41/s42/s37 JSONs untouched; compile clean at 15 pp, 0 errors (5 overfull
hboxes, cosmetic). The one open observation from round 2 (EXPERIMENT_PLAN.md's unfulfilled B3
promises removed rather than marked not-delivered) remains a no-action observation — the
deviation is documented in the audit.

## Score: 8.5 / 10

## Remaining weaknesses

No defects found this round. The remaining items are the paper's own disclosed limitations,
not defects, and each already carries its caveat in the text:

1. All membership-level arms are simulated or semi-synthetic (real 1000G is used with synthetic
   splits; the biobank cohort is mosaic-copying) — disclosed in Threat Model/Discussion.
2. The missingness channel is constructed-by-design (ρ = member-heavy batch fraction), with the
   ρ=0 null exactly at chance (0.502) as the control — disclosed in §scale.
3. s37/s41/s42 use 3 seeds/split per cell — tight stds on the headline deltas; a reviewer may
   ask for more seeds on the s37 relatedness control. Optional, not blocking.
4. ε_lb is framed as ordering information (valid flags disclosed) — correctly no longer
   presented as a calibrated absolute guarantee.

Minimum fixes: none required.

## Verdict

**YES — ready for submission** as a measurement paper. Every claim-vs-artifact mismatch found
across three review rounds is fixed and re-verified in the artifacts; the headline round-11
story (built sole-egress binding, bundle composition with a falsification control, multi-chrom
cM IBD bound with a planted control, ρ=0 null at chance, biobank-scale cohort with passing
controls) matches its JSONs to the quoted precision, and the paper's scoping statements are
honest throughout. Recommend: treat round 4 as a pure regression check; no further edits
demanded.
