# Cross-model review bundle — Audit-Governed Multi-Site GWAS (experiment results)

You are a hostile but fair reviewer (PoPETs/USENIX/CCS audience). The executor
ran CPU-only experiments for a paper on **auditable reveal governance for
multi-site GWAS** and needs a verdict: are these results sufficient to write a
submission-worthy paper, and if not, what is the single most important gap?

## What was run (all simulated, CPU-only, seeded; scripts in code/, raw in results/raw/)
- Two simulators: parametric Balding-Nichols multi-population panel; msprime
  3-population recombination-realistic demography.
- Attack families: (i) aggregate allele-frequency likelihood-ratio (Homer/
  Sankararaman); (ii) Beacon likelihood-ratio (Shringarpure/Bustamante);
  (iii) ancestry-only control attack (PCA-distance, confound upper bound).
- Controls: study members vs non-members drawn from the SAME population
  (ancestry- and frequency-base matched); attacker reference disjoint from study.

## Results (full detail in refine-logs/EXPERIMENT_RESULTS.md)
- S0 sanity: uncorrected λ_GC 6.35 → PC-corrected 0.979; attack AUROC 0.923;
  ancestry control 0.499; noise 0.483 → meter valid.
- S1 frontier (5 seeds): AUROC by release tier — F2 top-200 0.678; F4 all-SNP
  AF-1dp 0.562; F3 all-SNP AF-2dp 0.911; F5 all-SNP full-AF 0.916; F5-per-site
  0.994. Ancestry-only control ≈0.49 everywhere. Per-ancestor AUROC uniform.
  Site-size sweep TPR@FPR≤1%: 0.89 (n=50) → 0.68 (100) → 0.49 (200) → 0.43
  (300) → 0.27 (600/900).
- S2 intermediate reveals (5-fold CV logistic meta-attack, 3 seeds): final-only
  AUROC 0.912; +PCA-spectrum 0.949 (**Δ +0.037**); +λ_GC 0.000; +condition bit
  0.000; +missingness +0.001.
- S3 sanitized-diagnostics paradox (2 scripted analyst families × 2 seeds):
  probe reveals without sanitized diagnostics 3.0/4.0; with sanitized 0.0/1.0
  → 100%/75% reduction at EQUAL analysis quality (quality=5 both arms).
- S4 msprime cross-check (3 seeds): full-AF AUROC 0.716, AF-1dp 0.627; same
  direction as the parametric simulator.

## Proposed claims
C1 (credibility): individual-level membership leakage from released AF survives
ancestry-matched, panel-disjoint controls; ancestry-only control is at chance.
C2 (spine): intermediate PCA-spectrum reveals add measurable marginal leakage
(+0.037 AUROC) beyond the final summary release; λ_GC/condition-bit are
study-level and add ~0 individual leakage (reported honestly).
C3 (bet): pre-providing sanitized diagnostics reduces cumulative exposure by
cutting probe reveals (100%/75%), at equal quality.

## Known weaknesses (state your view on each)
- All simulated; no real cohort or real published summary-stat file.
- Beacon family ≈ chance at tested settings (structurally weak, not a true null).
- F1 (significance-only) tier degenerate (1 visible SNP).
- S3 analyst is scripted, not empirical.
- No MPC / enforcement prototype built — the "reveal governance" mechanism is
  framed but not implemented in this run.

## Your task
1. Verdict: are C1–C3 supported by this evidence? Mark each SUPPORTED /
   PARTIAL / NOT SUPPORTED with one-line justification.
2. The single most important experiment missing before submission.
3. Is "sanitized-diagnostics paradox" a real finding or a model artifact?
4. Rank the three claims by how much a reviewer would value them.
5. Name the strongest rejection argument and the cheapest experiment that
   pre-empts it. Do not propose adding modules; rank and judge what exists.
