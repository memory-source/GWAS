1. **The DP-inadmissibility map of GWAS artifacts** — Sensitivity of λ_GC, PCA spectrum, missingness, kinship submatrices, per-SNP p-value vectors. Contradiction: unbounded-sensitivity artifacts mean no DP budget admits them at GWAS utility. LOW risk, days–weeks.
2. **Marginal leakage of intermediate QC diagnostics** — PCA spectrum/loadings alone support above-chance site MI; λ_GC alone leaks little. LOW–MEDIUM, weeks.
3. **Precision non-monotonicity: does coarser output ever leak more?** — MEDIUM, days–weeks.
4. **Do release relations track attack-risk growth?** — union-coverage/refinement predict AUROC growth; receiver-expand does not. MEDIUM, weeks.
5. **Does the summary-stat leakage frontier survive matched controls?** — residual leakage concentrates in low-MAF tail and effect-size vectors. LOW, weeks.
6. **Manifest grammar coverage** — ≤5 constructs cover >90% of benign reveals. LOW–MEDIUM, weeks.
7. **Bypass-injection bake-off** — enforced egress vs logging-only provenance vs per-site API egress. MEDIUM–HIGH, months (1–2 if scoped).
8. **Slicing attack** — build forbidden matrix from authorized scalar reveals; sharp granularity threshold. MEDIUM, weeks.
9. **Statistical cost of pre-commitment** — MEDIUM, weeks.
10. **Where composition additivity fails for genotype attacks** — HIGH, weeks–months.
11. **Sanitized-diagnostics paradox** — releasing sanitized diagnostics reduces cumulative risk by removing probing incentive. HIGH, weeks.
12. **Receiver-dimension collapse** — most non-kinship artifacts become public. LOW, days–weeks.

(Full text verbatim in the harness tool result; condensed here for trace compactness.)
