# Deep-Review Log — Selection-Input Spectral Tomography (2026-10-16)

Direction: `research-plan/20261016_deep_search/DEEP_DIRECTIONS_ZH.md`
Backend: codex (gpt-6-astra, read-only, cross-family vs the deepseek-flash executor)
Thread: `01a14176-aaa4-7bf0-bf4c-9065eda1a1fe`
New artifacts: `refine-logs/EXPERIMENT_PLAN_DEEP_20261016.md`,
`code/s61_spectral_tomography.py`, `results/raw/s61_spectral_tomography_chr22.json`,
`refine-logs/EXPERIMENT_RESULTS_DEEP_20261016.md`, `EXPERIMENT_AUDIT_DEEP.md`/`.json`.

## Phase 1 — codex deep-review loop (converged, verdict ready)

### Round 1 — Score 4/10, Verdict not ready
Critical weaknesses (ranked):
1. The AF baseline was not given the same allowed information as the spectral
   attack (only the reference estimated the background frequency). → Add pooled
   background (reference+partner+shadow) and a conditioned single-release spectral
   baseline.
2. Scatter recovery is numerically fragile at the fixed `1e-9` cutoff (min singular
   value ~6e-12) and breaks under float32/6dp export. → Add an export-precision
   probe and limit the scatter claim to ideal precision.
3. The real workflow interface (plaintext loadings/eigenvalues + center/scale) is
   asserted, not verified against an actual pipeline.
4. The adaptive strategy shows no advantage; the coverage curves are partly
   tautological (visible coverage necessarily rises with more outputs).
5. Scientific cost is not measured (the GWAS row recomputes the honest PCA, not the
   per-release joint PCA).
6. Scope: one chromosome prefix, n=120/360 honest rows, no cross-cohort check.

### Round 2 — Score 5/10, Verdict not ready
Fixed: pooled/conditioned baselines added (information matched), second
within-population arm (EAS), export-precision probe, changed-count semantics.
New critical weakness: `spec_cond` fixed its orientation AFTER removing the
calibration members, so its calibration half was uninformative and the sign
flipped randomly per split; corrected `spec_cond` is EUR 0.985 / 3pop 0.757 /
EAS 0.987, lowering the EUR increment from +0.039 to +0.015.
Also flagged: the probe quantized the attacker's own `K` (should only quantize
released V/Λ), and the precision table showed seed0 only.

### Round 3 — Score 5/10, Verdict not ready
Fixed: `spec_cond` orientation fixed on the unconditioned release (matching the
conditioned-AF path); the precision probe quantizes only released quantities and
uses a precision-matched strong passive baseline; the report aggregates all five
splits and includes EAS in the GWAS range. Codex confirmed the corrected q8 gaps
and that the report matches the JSON.

### Rounds 4–8 — Score 5 → 6 → 6.5
Added per-release scientific cost (`sci_cost`), chr1 replication, real-software
export (PLINK 2.0), nested honest-sample-size ladder. Each new artifact was
independently recomputed and the report/JSON reconciled. Codex closed the
implementation and scope items in turn and raised the score to 6.5 ("almost").

### Rounds 9–10 — Score 7/10, Verdict ready (converged)
Round 9 closed the export-interface item (verdict "almost"). Round 10 closed three
closing items: the s63 `K`-shadowing artifact (408 MB JSON → 7.8 KB, `"k": 10`),
baseline consistency (s63 now reports conditioned AF), and the report wording. Codex
returned **verdict ready**: for a bounded measurement paper under the declared
information conditions, it is submittable — no further large experiment is required
as a submission gate.

### Status: CONVERGED — verdict YES (ready) at 7/10

Codex confirmed the round-10 closing items and changed the verdict from "almost" to
"ready": for a bounded measurement paper under the declared information conditions
and the real PLINK export precision, it is submittable. Score progression across the
ten rounds: 4 → 5 → 5 → 6 → 6 → 6.5 → 6.5 → 6.5 → 7 → 7. No further large
experiment is required as a submission gate.

## Phase 2 — pi-review (glm-flash, cross-family) loop

A fresh `/auto-review-loop` was run with the local pi agent as an independent second
cross-family reviewer. Thread `54e81425ae86484ba6c4e751d772fe15`.

### Round 1 — Score 6.5/10, Verdict almost
Seven weaknesses (W1–W7): W1 effects reported only as ceiling-censored AUC; W2 the
capability question at a single design point (q ≤ 8); W3 no-release control absent;
W4 adaptive framed as machinery; W5 effect heterogeneity under-narrated; W6 interface
generality asserted from one package; W7 nits (unused `Zh`, bootstrap seed, HWE
wording, span-restricted scatter). The reviewer independently ran an identity
ablation, a no-release control, and a ceiling analysis, and could not break the core.

Actions: added TPR@5%FPR columns to every table (s61/s62/s63); added the no-release
control (`cal_pca_k`/`cal_pca_qk`/`random_k`, |z| ≤ 5-clipped to remove a
preprocessing artifact); added the q-saturation cell (s64, q = 16 on three chr22
cohorts); reframed adaptive as a stronger-threat negative result; declared the
interface scope; removed the unused `Zh` param and documented the bootstrap seed and
the span-restricted scatter claim.

### Round 2 — Score 7/10, Verdict almost
W1/W2/W3 verified genuine (the reviewer re-ran the pipeline: the current code
reproduces the JSON to ≤1.4e-8 with bit-identical swap sequences; s64 matches
exactly). New findings: R1 the EAS random-span control 0.579 was a preprocessing
artifact (global-MAF panel + monomorphic markers), corrected to 0.463 with the q8 gap
surviving; R2 the W4/W6 sentences were claimed but not written (now written); R3 the
adaptive-CI list was incomplete (now all three cells); R4 s64 was unintegrated and
"engulfed" beyond the censored cohort (now integrated with the 3pop q16 headroom
cell); R5 TPR columns missing from some tables; R6 the fragile EAS-q4 claim.

Actions: clipped the control rows; corrected the EAS control to 0.463; ran s64 on
3pop (q8→q16 +0.049/+0.048, TPR 0.580→0.775); wrote the stronger-threat and
interface-scope sentences; added TPR to all tables; corrected the adaptive-CI list
and the effect-heterogeneity wording.

## Phase 3 — codex round-2 integrity re-audit

Thread `01a14455-c02b-7c93-b103-24ed8aace24a`, verdict FAIL → WARN → **PASS** within
the declared scope. It re-verified all new artifacts (105 s65 fields, 83 s66 fields,
322 table values) and closed: reporting contradictions (chr1 q4 5/5 → 3/5; adaptive
max cell; EAS-q4 list; clipping disclosure), a missing s64 tracker row, a
reviewer-only ablation (now artifact s65), an S66 coordinate-system bug (fixed), an
overstated identity-necessity claim (scoped), an overstated "all at chance" (now the
descriptive 0.455–0.509 bound), a mislabeled s65 field, and two rounding slips.

## Code-honesty note
The executor drove phase 1 and phase 3 with `mcp__codex__codex` (gpt-6-astra) and
phase 2 with `mcp__pi-review__pi` (glm-flash) — both cross-family backends, read-only
against the repo. Every numeric claim was independently recomputed by a reviewer; the
pi reviewer re-ran the pipeline (≤1.4e-8, bit-identical swaps), the codex auditor
replayed 105 s65 + 83 s66 fields and 322 table values.

## What was fixed (cumulative)
- Information-matched AF baselines (`af_ref`, `af_pooled`, `af_cond`).
- Conditioned single-release spectral baseline (`spec_cond`), orientation fixed on
  the unconditioned release.
- Export-precision probe over released V/Λ (float32, 6 dp, 4 dp) with a
  precision-matched strong passive baseline.
- `gap_vs_strongest_passive` with bootstrap CIs at every checkpoint.
- TPR@5%FPR columns (s61/s62/s63) to counter AUC ceiling-censoring.
- No-release control (`cal_pca_k`/`cal_pca_qk`/`random_k`), |z| ≤ 5-clipped, plus the
  s66 three-preprocessing provenance check.
- q-saturation cell (s64) to q = 16 on three chr22 cohorts.
- Identity ablation made an artifact (s65).
- Second within-population arm (EAS); second chromosome (chr1); nested
  honest-sample-size ladder (s62); real PLINK 2.0 export read-back (s63).
- Honest-only GWAS-power fit; per-release `sci_cost`; membership score free of a
  recovered-rank divisor; permutation null on the evaluated half; adaptive ranking
  free of default `pinv`.
- Artifact fixes: chromosome-scoped output filenames; the s63 `K` shadowing that
  wrote a 408 MB JSON; the s62 role truncation replaced by an explicit error plus a
  nested design.

## Terminal headline (bounded)
Under the declared information conditions, accumulating several final outputs of an
AF-preserving collaborative PCA beats the strongest passive baseline (composition-
balanced conditioned AF and the conditioned single release) by a small membership
increment (+0.014 [+0.007,+0.021], 5/5) that survives a real PLINK 2.0 export at its
printed precision; the practical reading is TPR@5%FPR 1.000 (EUR/EAS) and 0.580
(3pop) at q8; the full scatter recovery remains a float64-only witness.

## Method Description
The direction measures ACTIVE leakage from a *sequence* of final outputs, not a
single fixed release. A cooperating analysis partner holds an "own" genotype panel
and varies its own input across re-fits of a joint PCA with the honest study. Each
re-fit returns a truncated final output — the top-k joint eigenvectors `V_t` and
eigenvalues `Λ_t` — which the partner (who also knows its own scatter `K_t = Z_P(t)^T Z_P(t)`)
converts into `Y_t = V_t Λ_t − K_t V_t`. The central identity `Y_t = C_H V_t` is exact:
each output is a linear measurement of the honest sample scatter `C_H` along the
released directions. Stacking q outputs gives `C_H Q = Y` over `Q = [V_1..V_q]`, so the
range of `Y` recovers the span-restricted operator `C_H P_Q` and the projection norm
onto that range is a study-inclusion score. Input variation is an AF-preserving block
exchange between two of the partner's own rows, which leaves every per-SNP histogram,
allele frequency, variance and HWE/missing count unchanged while moving individual
structure. Three component studies make the mechanism verifiable, and the adversarial
ablation in s65 (raw-eigenvector stacking reproduces nothing) plus the no-release
control in s66 (known calibration-member genotypes at chance under three preprocessing
choices) isolate the effect to the tomography identity and the release, not to the
partner's own data. Data flow: real 1000G chr22/chr1 prefixes → disjoint 5-role split
(honest/partner/non-member/reference/shadow) → joint-PCA re-fits → stacked recovered
operator → held-out-half-oriented membership AUROC and TPR@5%FPR vs information-matched
passive baselines (conditioned AF, conditioned single release) → per-split JSON →
results tables.



