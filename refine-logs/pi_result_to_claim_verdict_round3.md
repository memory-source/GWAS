# Result-to-Claim Evaluation — Round 11 (C1–C5)

**Date**: 2026-10-07. **Method**: every cited number re-verified by direct JSON parsing of
`results/raw/s3[789]_*.json`, `s4[012]_*.json`, `s23_governed.json`, `s37_realsynth_scale.json`
(stronger than `grep -F`); paper claim paragraphs read in `paper/main.tex`; claim wording taken from
`refine-logs/EXPERIMENT_PLAN.md` Claim Map. No results were trusted from summaries alone.

## 1. claim_supported

| Claim | Verdict | One-line basis |
|---|---|---|
| C1 per-individual re-expression (A\*, ρ, f-DP ε_lb; control at ~0) | **yes** | All minimum evidence delivered in s38; ordering preserved; controls at 0.03/0.04; two disclosed caveats (below) |
| C2 real-data PC<AF is not a size-cap artifact | **yes** | s37: AF leads at every cell nsp 150→12k, gap plateaus at −0.03..−0.04 while in-sample gap →0; controls pass |
| C3 sole-egress binding holds | **yes** (scoped) | s39: bypass 3/3 → blocked 3/3; declared flows; chain verifies; transported tier byte-identical to s23 seeds 0–2 |
| C4 bundle > max single | **yes** | s40: +0.015 in 5/5 seeds; structural-zero control pinned (Δ −0.0003); modest, config-specific magnitude |
| C5 second channels bounded on realistic multi-chr data | **yes** | s41: IBD ≤ chance ≤ AF on 3 chromosomes with validated planted control; s42: strictly monotone in ρ, ρ=0 at chance |

## 2. what_results_support

- **C1** (s38, 5 seeds): A* pc_loadings 0.610 > af_full 0.546 > af_1dp 0.187 matches the AUROC
  ordering 0.882 > 0.845 > 0.622 (af_sitewise 0.805 tops both). pc_refbasis control A* 0.033 /
  ε_lb(0.01) 0.19; pc_eigonly A* 0.043 / 0.24 — "control at ~0" met. A*−(2AUC−1) = −0.056..−0.153 on
  all five leaking channels, +0.029/+0.043 on the controls: the metric is demonstrably not 2AUC−1 in
  disguise, and the paper reports the actual direction ("below it for the leaking channels").
- **C2** (s37, 3 seeds/cell): PC−AF gap −0.096/−0.062/−0.066/−0.030/−0.039/−0.034 across
  nsp 150/500/2000/8000/10000/12000 (n_per_group 450→4000; largest cell 36,000 study genomes): AF
  leads everywhere and the gap plateaus rather than closing. In-sample gap 0.199→0.081→0.025→0.0002.
  Controls: LD decay matches donor (incl. NaN-pattern match at nsp12000), relatedness within ~10% at
  5/6 cells (nsp150 +10.6%), AF spectrum within ~0.001, eig10 0.211 vs 0.218. Real haplotype
  structure from the 1000G phased library; combined arm ≡ af_only out-of-sample.
- **C3** (s39 + s23): library_leak 3/3 → natural_egress_blocked 3/3 (leak file never created,
  blocked `open()` recorded), sanctioned-undeclared blocked in both arms, declared flow 3/3, chain
  verified 3/3. Governed per-seed AUCs 0.62576/0.61781/0.62817 are byte-identical to s23 seeds 0–2
  (agg 0.6239 vs s23 0.6217) — "lands at the safe tier" verified exactly.
- **C4** (s40, 5 seeds): bundle 0.9619 > max_single 0.9471 (argmax per-site AF), excess +0.0149
  (per-seed +0.011..+0.018, 5/5 positive). af_plus_noise 0.8389 vs af_only_bundle 0.8392 (Δ −0.0003):
  feature count alone adds nothing. LOO: per-site AF drop 0.076; pooled AF + call-rate supply the
  excess.
- **C5** (s41, 3 splits; s42, 3 seeds/cell): real phased IBD 0.456/0.479/0.424 (chr1/2/22) ≤ chance,
  far below AF 0.866/0.864/0.838; planted 10-cM control detected on every chromosome
  (0.715/0.729/0.702; planted segment 9.98 cM) — the null is an honest bound, not a broken detector.
  Missingness strictly monotone in ρ (0.502/0.641/0.762/0.891/1.000), ρ=0 exactly at chance, bonus
  sweep 0.502→0.580→0.848→1.0, positive controls (superpop keying, no-holdout) at 1.0.

## 3. what_results_dont_support

- **C1**: the plan's binormal-direction success criterion (A*−(2AUC−1)>0) **failed** — the gap went
  negative on every leaking channel (honestly re-reported; still evidences non-binormality). Strict
  empirical f-DP bounds are **infeasible (valid=False) on all four high-AUROC tiers at every tested
  δ** (zero observed FNs); the quoted 4.997/4.12 rest entirely on the 1/900 resolution floor. No
  per-individual numbers exist on real 1000G or the s37 cohort (sim-only scope).
- **C2**: no cell revives PC>AF — "reproducing the crossover" holds only as the AF-first side
  persisting at all n; the crossover itself is located by the Fst sweep (measured ≈0.06 vs real
  ≈0.10), not by s37. Largest delivered cell is nsp=12,000 / 36k study genomes (~14× the 1000G arm),
  **not** the plan's nsp=33,000 biobank cell (cut at 24 GB GPU, deviation documented) and not
  UKB-scale; the out-of-sample gap narrows 3× from nsp150→nsp8000 before plateauing, so
  extrapolation to 500k is not directly measured.
- **C3**: nothing in the results contradicts the claim, but support is confined to the stated threat
  model (careless Python author, no native code, no root): PEP-578 hooks are advisory against native
  code, no adversarial bypass was attempted, single platform, demo not proof.
- **C4**: magnitude is regime-dependent: +1.5 AUC points on one simulated manifest, 5 seeds;
  max_single is oracle-selected on full data (makes the excess conservative, but the number is not a
  general constant — at s38's uncoarsened sitewise-AF point the single channel reaches 0.964,
  comparable to the bundle).
- **C5**: 3 of 22 autosomes (disclosed); 50k uniformly-strided markers, not full-density (disclosed);
  the missingness channel is constructed by design (ρ is the designed batch-composition parameter,
  chr22 only) — ρ=0 at chance shows the attack reads composition, not that an observed real-world
  batch leak was measured; base-missingness and sigma_batch sweeps are ceiling-saturated (no graded
  information there; disclosed).

## 4. missing_evidence

1. C2: one cell at nsp ≥ 33k (memory-lean int8 or CPU-sharded) to extend the −0.03 plateau beyond 12k.
2. C1: per-individual A*/ρ/ε at one real-1000G or s37 operating point; a confidence-corrected FNR
   bound (rule-of-three / Clopper–Pearson) so ε_lb is quotable beyond ordering information.
3. C4: bundle vs an uncoarsened per-site AF single channel; one alternative manifest/dataset.
4. C5: the remaining 19 autosomes for IBD; a real observed batch-effect event would upgrade
   "realistic" from skeleton-level to observed.
5. C3: nothing essential; optionally a documented native-code write attempt to characterize the
   threat-model edge.

## 5. suggested_claim_revision

- **C1**: keep. Fix the paper clause "so it is resolvable for δ≳0.01" (s38 shows valid=False even at
  δ=0.05/0.01 on af_full/af_2dp/af_sitewise/pc_loadings): state that on the four high-AUROC tiers the
  strict empirical trade-off is infeasible at all tested δ and the values are floor-driven ordering
  information; add the per-δ validity flags to the paper (currently only in EXPERIMENT_RESULTS B1).
- **C2**: keep the conclusion, fix the magnitude language: replace "With real haplotype structure and
  biobank-scale $n$" (main.tex ~line 773) with the delivered figure ("up to n_sp=12{,}000 — 36{,}
  000 study genomes, ≈14× the 1000G arm"); note the gap plateaus at −0.03..−0.04 rather than closing.
  Limitations already says "n_sp up to 12000" — align the scale-section sentence with it.
- **C3**: no revision; scoping is already honest ("ordinary Python, not root, no native code").
- **C4**: optionally add "for this manifest" and note max_single is oracle-selected (conservative for
  the claim).
- **C5**: rename the section header "bounded at physical length" → "bounded in genetic (cM) length"
  (the body already says centiMorgans via the IMPUTE2 map); the body's "a real function of batch
  composition rather than a construction artefact" can be softened to "driven by the batch-composition
  parameter rather than a skeleton artefact".

## 6. next_experiments_needed

Ranked by claim-hardening per cost:
1. **s37 extension** (C2): one nsp≈33k cell — directly tests the plateau beyond 12k.
2. **ε_lb recompute** (C1): confidence-corrected FNR floor as the headline bound + paper disclosure
   of validity flags.
3. **Per-individual metric on real/s37 data** (C1): one operating point.
4. **IBD on remaining autosomes** (C5): pipeline exists; parse cost only.
5. **s40 robustness** (C4): uncoarsened per-site AF single channel + second manifest.
6. Optional **C3 boundary test**: native-code write attempt, documented as out-of-model.

## 7. confidence

**high** — all 20 cited values reproduced by direct JSON parsing (not just `grep -F`); per-seed data
internally consistent (s39↔s23 byte-identical on matching seeds); paper wording checked against
artifacts; the three residual issues (ε_lb resolvability clause, "biobank-scale n", "physical length"
header) are text-level, disclosed-or-minor, and do not change any verdict. Residual uncertainty: C2
rests on a plateau extrapolated from 12k, and C3 is a scoped demonstration, not a proof.
