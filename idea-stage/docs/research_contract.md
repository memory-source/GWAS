# Research Contract — Audit-Governed Multi-Site GWAS

Selected idea: **T5 + T2 + T11** (matched-control leakage frontier; marginal intermediate-reveal leakage; sanitized-diagnostics paradox), from `idea-stage/IDEA_REPORT.md`.

## Primary claims (each must be traceable to a result file)
- **C1 (measurement credibility / T5):** On the F5 (full-summary-statistics) endpoint, at least one attack family achieves TPR@FPR≤1% significantly above chance, *and* an ancestry-only control attack is demonstrably weaker than the real attack (ordered by AUROC / TPR@low-FPR). Under ancestry-matched, panel-disjoint conditions the residual individual signal either (a) survives — locating leakage in the low-MAF / effect-size objects — or (b) vanishes (documented negative).
- **C2 (spine / T2, proposal C3):** Adding intermediate diagnostic reveals (λ_GC; PCA spectrum + PC–phenotype correlations; missingness/QC tables; condition bit) on top of a fixed final-release tier produces a direction-consistent, measurable marginal leakage increment for at least one diagnostic type, at TPR@FPR≤1% under ancestor-matched, panel-disjoint controls, in at least one final-release tier. (Both signs valid.)
- **C3 (bet / T11):** Pre-providing sanitized versions of commonly-probed diagnostics reduces the number of probe reveals an analyst issues at equal analysis quality by ≥20%, consistent across two independent scripted policy families and ≥2 seeds. (Dropped as policy artifact otherwise.)

## Non-claims (explicitly out of scope — write into the paper)
- No data-independent privacy bound / no ε-budget. Leakage numbers are model-specific estimates.
- No claim of monotone precision⇒leakage relation; no claim that any single scalar measures genomic privacy (Wagner 2016).
- No end-to-end security theorem; C1 of the proposal is a trace-compliance property on generated artifacts under a stated TCB, not a DP composition theorem.
- No claim about adversarial authors who can execute arbitrary code (outside A2).

## Threat model
- A1 honest-but-curious site (handled by MPC; not a contribution).
- A2 constrained workflow author (missing/over-precise reveals, wrong receiver, loose condition).
- A3 adaptive analyst (cross-run accumulation; detection only, no formal composition bound).

## Kill criteria (pre-registered)
- If S1's attack endpoint fails (no family beats chance on F5, or confound control is not weaker) ⇒ ALL leakage conclusions blocked; fix attacks first (proposal P0-3 fail branch).
- If S2 is negative ⇒ reframe to "release-strategy empirical frontier"; keep intermediate reveals as documented negative.
- If S3 fails the drop rule ⇒ drop (recorded).

## Evidence files
- Results → `results/raw/`, figures → `paper/figures/`, tables → `results/tables/`.
- Every claim above cites a specific JSON/CSV under `results/`.
