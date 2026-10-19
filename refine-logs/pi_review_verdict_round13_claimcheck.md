# Round-13 Result-to-Claim Evaluation (C6/C7/C8) — claim check

**Reviewer**: pi (independent artifact-level re-verification; all numbers below were re-read from the JSONs, code, and current `paper/main.tex` in this session).
**Date**: 2026-10-07
**Claim set**: C6 (joint size×structure governance + real-genome reproduction of the PC-dominant small-(n×M) regime), C7 (within-pop replication ≥3 chromosomes / ≥2 populations), C8 (in-sample loadings mechanism; held-out-half transfers at chance).
**Artifacts read**: `results/raw/s48_real_structure.json`, `s48_scalematched.json`, `s49_real_mechanism.json`, `s50_real_structure_multichrom.json`, `s21_sensitivity.json`, `s16_multichrom.json`; `code/s48_real_structure.py`, `s49_real_mechanism.py`, `s50_real_multichrom.py`, `code/common.py`; `paper/main.tex` (abstract, contributions, §scale, §boundary, Fig.7 caption); `refine-logs/EXPERIMENT_RESULTS.md` B2/B3/B4.

## Verdict: **claim_supported = YES** (C6 yes-with-scope, C7 yes, C8 yes). Confidence: **high**.

Every claimed number was re-verified against the artifacts; all match. The residual issues found are precision/disclosure items, none direction-threatening.

## 1. Evidence verification (all ✓)

- **s48 chr22 ladder** (5 split seeds, means): 12 single-pop cohorts gap +0.0502…+0.1339 (n_study 28–36) ✓; 5 two-pop mixes +0.1115…+0.1719 (n 52–66) ✓; 3pop −0.0053 @ n=450 (Fst 0.0874) ✓; 5pop −0.0296 @ n=715 (Fst 0.0747) ✓. Quoted AF/PC pairs (GBR 0.896/0.998, CHB 0.863/0.996, TSI 0.861/0.995, CDX 0.884/1.000) exact ✓.
- **s48 scale-matched** (`cohort_cell(per=20)`): GBR +0.091, CHB +0.0765 @ n=20; 3pop +0.1607 @ n=60, Fst 0.0997 ✓; 5pop +0.1112 @ n=100, Fst 0.0885 ✓. Controls refb 0.5195–0.5625, eigo 0.5182–0.5965 ✓.
- **s49 mechanism** (chr1/chr2 × GBR/CHB/YRI, 6 cells, 5 seeds): top3 0.881–0.959 ✓, full 0.9943–1.0 ✓, eig_only 0.5447–0.5964 ✓, refbasis 0.5516–0.5951 ✓, oos_transfer 0.5512–0.5809 ✓, samehalf 0.9981–1.0 ✓, n_study 30–36, k=10 ✓.
- **s50 chr1/chr2**: within-pop +0.0511…+0.0858 all positive ✓; across-pop −0.0391…+0.0231 with chr2 3pop +0.0231 exception ✓ (Fst 0.0855, n=450).
- **Code audit**: disjoint 3-way same-population splits (negatives never from the attacker's own reference panel — the round-1 flaw stays fixed); PCA fit on study only; mean centering from study; `auroc` tie-safe; 5 seeds; s50 is a faithful re-run of s48's `run_chrom`. Methodology sound.

## 2. What the results support

- **C7 — fully.** Within-pop PC>AF on 3 chromosomes (chr22 × 12 populations, chr1 × 3, chr2 × 3 = 18 cohort-cells, all positive, 5 seeds each) plus n=20 scale-matched GBR/CHB. "≥3 chromosomes / ≥2 populations" is exceeded.
- **C8 — fully.** The control quartet cleanly isolates the mechanism: study-fitted loadings ≈1.0; random-basis (eig-only) 0.545–0.596; reference-panel basis 0.552–0.595; held-out-half basis transferred to the other half 0.551–0.581 (chance); same-half refit 0.998–1.0. In-sample-fitted directions are necessary and sufficient — matching the simulator mechanism (ref-fit collapse, oos at chance).
- **C6 — yes, with an explicit division of labor.** (i) Structure alone is refuted on real data: at matched small n the PC channel dominates even at Fst≈0.09–0.10 (scale-matched +0.161/+0.111). (ii) Size gates the switch: at Fst≈0.075–0.10 the gap goes ≈0 at n=450 (−0.005…+0.023) and clearly negative at n=715 (−0.030…−0.039, all 3 chromosomes). (iii) The structure axis at fixed size (monotone decay, crossing rising as n shrinks) is carried by the simulator (s21: measured-Fst crossing brackets (0.027–0.052) @n4800, (0.052–0.077) @n2400, (0.103–0.151) @n900, 15 seeds — which is *stronger* joint-law evidence than the text quotes). (iv) The small-(n×M) PC-dominant regime is reproduced on real genomes within single populations.
- W13-1…W13-5 fixes verified as landed (true s50 range + chr2 exception; "monotone" confined to the size-controlled sweep; sizes + scale-matched rung reported; k/n≈1/3 disclosed; §boundary sentence parses; results-doc B2 now chr1/chr2).

## 3. What the results do NOT support

- **No real low-structure large-n cell.** Whether AF can beat PC at Fst≈0 given size (AND-gate vs shifted-threshold form of the joint law) is undetermined on real data; the ladder's large-n rungs all have structure ≥0.075.
- **No uncertainty on real gaps.** JSONs store only 5-seed means (no per-split values, CIs, or permutation p-values). The 3pop cells are quoted as "one point beyond the crossing", but the gap is ≈0 with chromosome-level sign instability (−0.005 chr22/chr1 vs +0.023 chr2); "beyond" is unsupported — "at the crossing (≈0)" is the honest statement.
- **M-axis on real data is one thinning contrast**: chr2 flips (−0.063 @ M=4664 → +0.023 @ M=3201) but chr1 only shrinks (−0.09 @ M=5257 → −0.0055 @ M=3384); all s48–s50 runs hold M≈3200–3391. The real n-axis is well covered; the real M-axis is not.
- **Gap magnitudes are not comparable across ladder rungs** (k/n runs 1/3 → 1/70 with n; only signs/orderings are meaningful). Disclosed for s49, not restated for s48's large-n cells.
- **PC-arm orientation uses test labels** (`_cal = max(a,1−a)` in `cohort_cell`/`cell`), while the AF arm is never flipped and the simulator protocol uses held-out-half orientation. This inflates null-level PC values (explains the 0.51–0.60 control band) and adds ~±0.016 noise at near-chance cells — the same order as the 3pop gaps. Headline cells (0.99+ vs 0.86–0.95) are unaffected.

## 4. Missing evidence

1. Per-split values / bootstrap CIs / permutation p-values for real cells (`stratified_permutation_auroc` already exists in `common.py` — unused in s48–s50).
2. Within-pop n≈54–100 rung (YRI n=108 admits a 2-way split; currently capped at 36 by the 3-way split).
3. Intermediate-n structured rungs (3pop at per∈{50,100} → n=150/300) to trace the real crossing in n at fixed composition.
4. Held-out-half orientation re-score of s48/s50 (or disclosure of the train-calibrated orientation for real runs).
5. An s21 rung at nominal Fst 0.06 to tighten the sim bracket actually quoted.

## 5. Suggested claim revision (all minor; text-level)

- Abstract/contribution/Fig.7: replace "crosses at measured ≈0.04–0.06" with the per-size measured brackets (0.03–0.05 / 0.05–0.08 / 0.10–0.15) — this *strengthens* the joint law (crossing rises as n shrinks) and removes an interpolation.
- §scale: "negative at Fst≥0.075" → 5pop measured Fst is 0.0747; use "≳0.07" or the actual range.
- §scale: "one point beyond the crossing" → "at the crossing (gap ≈0 ± 0.02 across chromosomes)".
- §scale: "the unconfounded 3/5-pop cohorts" → they are confounded (both high); the *joint design* (ladder + scale-matched) does the unconfounding.
- §scale: "flips sign between M=4664 and M=3384" → true for chr2 only; chr1 shrinks without flipping. Reword to "moves toward PC as M falls (chr2 flips; chr1 approaches zero)".
- Surface the s50 chr1:YRI eigonly 0.618 (or state the s50 control range) — currently no quoted range covers it.
- Results-doc B2 "eigenvalue-only 0.53–0.60" → actual min 0.545 ("0.54–0.60").

## 6. Next experiments

1. **s51 within-pop large-n**: YRI/CHB-scale populations, 2-way split (study ≈54, ref/non ≈27), 5 seeds — fills (Fst≈0, n≈54); decides whether the gap decays toward 0 with size at zero structure.
2. **s52 fixed-composition n-ladder**: 3-super-pop cohorts at per∈{20,50,100,150} (n=60/150/300/450), 5 seeds — traces the real crossing in n at fixed Fst≈0.09; predicts a sign flip between 100 and 450.
3. **Re-score s48/s50** with held-out-half orientation + per-split logging + permutation p-values.
4. Optional: s21 rung at Fst 0.06; M-contrast run at matched (n, Fst) to give the M-axis one clean real cell.

## 7. Confidence

**High.** Artifact-level verification passed end-to-end; replication breadth (18 within-pop cells over 3 chromosomes/12 populations) and the mechanism quartet are decisive for C7/C8; C6's joint-governance framing is supported with the scope conditions in §3, and the simulator carries the crossing surface the real data cannot reach at biobank n.
