# Round-1 review verdict (pi-review / glm-flash)

The reviewer verified headline numbers against results/raw/*.json and read all
five scripts.

## Verdict
- **C1 (credibility) — PARTIAL.** Direction holds (ancestry control 0.492,
  per-ancestor uniform), but `attack_af()` uses the attacker's own reference
  panel as non-members → symmetric two-panel distinguishability (≈2× separation
  of the canonical one-sided test). Beacon family (0.511) is null *by
  construction* (common-SNP queries, genome-mean f → p0≈1).
- **C2 (spine) — PARTIAL.** +0.037 solid within setup (every +PCA seed exceeds
  every base seed; CV clean). Same negative-set bias applies; the zeros for
  λ_GC/condition-bit are trivial by construction (constant feature).
- **C3 (bet) — NOT SUPPORTED.** Probe counts are hardcoded constants;
  analysis_quality is arm-independent; cumulative_auroc is a max, not
  cumulative; noisy branch never fired; conviction rule guaranteed to pass.

## Answers
- Most important missing experiment: re-run S1/S2 with a THIRD disjoint
  non-member panel (pre-registered, cheap).
- Paradox: model artifact — salvageable only with a measured stopping rule.
- Ranking: C2 > C1 > C3.
- Strongest rejection: "the frontier is measured on a self-referential control";
  cheapest pre-empt: the disjoint-negative re-run.
- Hygiene: "monotonic" false at tail (0.267→0.274); bootstrap CI promised but
  absent; permutation decomposition function unused; S0 min_p_causal NaN.
