# Round-11 Review — Loop 2 (run_20261006_f68b0cd6)

**Reviewer**: pi (artifact-level verification; diffed `paper/main.tex` against the round-1
snapshot `main_20261006_232000.tex` — exactly six change hunks, all checked; re-read
`code/s38_privacy_metric.py`, `code/governance.py`, the re-run `s38_privacy_metric.json`, and
the three docs).
**Date**: 2026-10-07

---

## W1–W8 scorecard (round-1 weaknesses)

| # | Round-1 weakness | Status | Evidence |
|---|---|---|---|
| W1 | Paper: "A* exceeds 2AUC−1 for the leaking channels" (inverted vs JSON) | **FIXED** | Paper now: "A* differs from the binormal value 2AUC−1 (below it for the leaking channels, above it for the at-chance controls)" — matches `gap_vs_2aucm1` exactly (−0.056…−0.153 leaking, +0.029/+0.043 controls). |
| W2 | ε_lb: b1-only formula in text, clip-driven ≈24 values, n_neg≈300 | **FIXED** | Code: `eps_lower_bound(fpr,tpr,δ,n_pos,n_neg)`, FNR floored at 1/n_pos (docstring honest: "a single missed member is possible given the sample even when the empirical TPR is 1"); paper states the two-branch max and 1/n_neg=1/n_pos=1/900 with the ordering-information caveat; JSON re-run: pc_loadings 4.997→"≈5.0" ✓, af_1dp 1.057→"1.1" ✓; results doc updated. Estimator verified legitimate (both branches ≡ the two branches of g_{ε,δ} under the empirical-FNR convention). |
| W3 | §replication "remains future work" vs §scale reporting it done | **FIXED** | Now: "is reported in Section~\ref{sec:scale}, where it is bounded on three chromosomes at physical length." |
| W4 | Limitations cross-ref sec:governance for the cohort arm | **FIXED** | Now sec:scale. |
| W5 | "base missingness shifts the absolute AUROC only weakly" over a saturated sweep | **FIXED** | Now: "within the tested range the base missingness level leaves the channel saturated at ρ=1" — accurate (sweep 0.9992–1.000). |
| W6 | "at full marker density … strided uniformly" self-contradiction | **PARTIAL** | Parenthetical now discloses "50k markers per chromosome, uniformly strided from the full-density parse" — but the main clause still says "at full marker density", which the parenthetical contradicts (50k uniform stride ≠ full density; parse capped at 200k). |
| W7 | Audit claimed design-doc refresh not applied; CV 0.32 garble | **MOSTLY FIXED** | `ROUND11_DESIGN.md` now states the delivered 6-cell grid with the historical note (nsp=33000/thin-dense "exceeded the 24 GB GPU at the Gram-eigenvector step") ✓. Residuals below (R2, R3, R5). |
| W8 | s38 dead import; governance.py docstring "remains argued" | **FIXED** | `anc_only_attack` gone from s38; docstring now "The sole-egress binding over this monitor is built in `code/enforce.py` …". |

## Residual weaknesses (new/remaining, ranked)

**R1 (minor — paper).** §governance: "versus $0.53$/$4.2$ for full-precision AF" — these are
**af_2dp's** values (JSON: 0.5318 / 4.190), not af_full's (0.5456 / 4.120 → 0.55/4.1). The
results doc has it right ("af_full .845/.546/4.1; af_2dp .840/.532/4.2"); the paper's label is
wrong. Reviewer note: this mismatch existed in round 1 (there it hid behind "0.53/24" where 24
matched both tiers) and my round-1 weakness list missed it — caught now because the new ε
values disambiguate the tier.
**Minimum fix**: "0.55/4.1 for full-precision AF" or relabel "for AF@2dp".

**R2 (minor — docs).** `EXPERIMENT_RESULTS.md` B5: "(nsp150 141 vs 128 markers, +11%; **donor**
CoV 186/141=1.32 vs **cohort** 50/159=0.32)" — both attributions garbled: 186/141 is the
nsp150 **cohort** std/mean (cohort CV 1.32; the donor at nsp150 is 48/128 = 0.38), and 50/159
is the **nsp10000** cohort's CV (a different cell). The round-1 error (cohort CV quoted as
0.32) is now arithmetically present but mis-attributed.
**Minimum fix**: "cohort CoV 186/141=1.32 vs donor 48/128=0.38 at nsp150; other cells ~0.3".

**R3 (minor — docs).** `EXPERIMENT_TRACKER.md` S38 row still carries the round-1 wrong-direction
shorthand "A*−(2AUC−1)>0 (non-binormal)" — true only for the controls; the paper and results
doc are now correct.
**Minimum fix**: "<0 on the leaking channels, >0 on the controls" (or drop the sign clause).

**R4 (minor — paper, = W6 remainder).** s41 main clause still says "at full marker density".
**Minimum fix**: "at uniform genome-wide density (50k markers per chromosome)".

**R5 (optional).** The re-run JSON's `eps_lb_0.01_valid` is **False** on the four leaking
channels (the strict empirical trade-off is incompatible with any finite ε at the zero-FN point
— exactly what the 1/n_pos floor truncates). The paper's "ordering information rather than an
absolute guarantee" caveat covers the spirit; a one-clause disclosure in the results doc would
make the flags self-explanatory.

**Observation (no action demanded).** `EXPERIMENT_PLAN.md`'s unfulfilled B3 promises (MLP
cross-check, anc_only_auc) were removed rather than marked not-delivered. The deviation is
documented in the audit, so this is acceptable for an internal planning doc; pre-registration
practice would keep the original text plus a deviation note. The plan's nsp=33000 line is fine
as historical intent given the design-doc note.

## Verified unchanged (regression)

- Paper diff vs round-1 snapshot: exactly the six intended hunks; every other number identical.
- s38 A*/AUC/gap values identical between the round-1 and round-2 JSONs (same seeds/attacks;
  only eps_lb changed, as expected — A* does not depend on the estimator). eps values moved
  exactly as the FNR-floor predicts (af_full 24.01→4.12, pc 23.93→5.00, af_sitewise
  26.29→5.78; controls unchanged since their TPR never reaches 1; af_1dp unchanged at 1.057).
- s39/s40/s41/s42/s37 JSONs untouched; compile clean at 15 pp (5 overfull hboxes, cosmetic).

## Score: 8.0 / 10 → 8.5 after R1–R4 (four one-line fixes)

## Verdict

**Almost — one line-edit pass away.** Five of eight round-1 weaknesses fully fixed, one
partially (W6), and the residuals are R1 (one mislabeled tier in one sentence), R2 (one garbled
doc line), R3 (one stale tracker shorthand), R4 (one phrase). No new experiments; no new
contradictions introduced by the edits. The substantive round-11 story is artifact-backed and
now correctly stated: the two-branch ε_lb is honestly framed as ordering information, the A*
direction matches the JSON, the IBD future-work contradiction is gone, and the design doc
records the delivered grid. Apply R1–R4, recompile, and round 3 should be a pure regression
check.
