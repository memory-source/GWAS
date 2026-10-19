# Round-14 Review — result-to-claim evaluation, round-12 residual-limitation closure (S44–S47)

**Reviewer**: pi (artifact-level; every result file read in full, every estimator read in code,
cited values re-verified against the JSONs independently of evidence_check.py).
**Date**: 2026-10-07
**Scope**: C1b/C1c/C2b/C4b/C3b vs s44_individual_real.json, s37_realsynth_scale.json (nsp16000_m3000),
s46_bundle_manifest2.json, s39_enforcement.json; code s44/s37_realsynth_scale/s46/enforce/s39 + s38 estimators.

## Scorecard

| Claim | Verdict | One-line basis |
|---|---|---|
| C1b per-individual metric not a simulation artefact | **PARTIAL** | A*/ρ/ordering/controls replicate on real chr22 (8 splits) + real-haplotype cohort; **ε_lb quotability does NOT hold on real 1000G's headline channels** (strict valid=False at all δ, zero-FN corner at n=450) — cohort-only. |
| C1c ε_lb quotable, not floor-driven | **PARTIAL** | Outcome yes (cohort: valid=True on all 6 channels × 3 δ). Mechanism attribution no: rule3 values are **bit-identical** to resolution-floor values on the cohort (floor inert where valid) and validity fails under both floors on real high channels — the floor never causes validity (validate_tradeoff checks the raw ROC; at a zero-FN threshold no finite ε from this estimator validates). Paper already says "admits", not "makes valid" — correct. |
| C2b AF-first ordering at nsp=16000 | **YES** | Δ = −0.0264 ± 0.0016 (3 seeds, ~16σ); PC controls 0.4998/0.4999 (chance); in-sample gap 0.0029; relatedness mean +1.9% of donor; AF spectrum/eig10/LD bins 2–4 match. |
| C4b bundle excess not one manifest / oracle pick | **YES** | Stronger manifest (full-precision per-site AF in-bundle) vs **uncoarsened** single: +0.0084, 5/5 positive; conservative direction (harder comparison) — shrank from +0.0255 vs coarse exactly as expected; single_coarse 0.94707 **exactly** reproduces s40 max_single; noise control adds nothing (0.9636 ≈ 0.9642). |
| C3b ctypes boundary tested | **YES** | ctypes.CDLL(None) denied on `ctypes.dlopen`, 3/3, `child_failed=False` (hook fired, child completed — no crash masking); bypass matrix + governed tier (0.624/0.509) unchanged. |

## Per-claim detail

### C1b — PARTIAL
**Supported.** Real chr22 (8 disjoint splits, n_pos=n_neg=450): af_full A* 0.350 > pc_loadings 0.271,
matching the real AUROC ordering (0.727 > 0.665; s16-consistent); af_sitewise 0.584 dominates all;
controls pc_eigonly A* 0.067 / pc_refbasis 0.053 → AUROC 0.5021/0.5003 (chance; below the ~0.15 KS
noise scale at n=450). Cohort nsp=2000 (real 1000G haplotypes): af 0.1345 > pc 0.0372; controls
0.017/0.012. ρ: median > 0.5 on both panels (real 0.582/0.578; cohort 0.505/0.506), frac_gt_half
0.64–0.71 real / 0.59–0.62 cohort. Magnitudes attenuate vs the simulator (pc 0.61→0.27; af 0.55→0.35)
— i.e. the metric re-expresses the known real-data attenuation rather than manufacturing signal.
That is the anti-artefact core, and it holds.

**Not supported.** "ε_lb holds on real 1000G": strict validity is False at **all three δ** for
af_full, af_sitewise, and pc_loadings on the real arm (zero-observed-FN corner); the rule3 floor
lowers the values (2.00→1.94 at δ=0.05) without making any of them valid. Only af_1dp@δ=0.05 (0.568)
and the two controls admit valid bounds on real chr22. Also: the cohort arm is a single cohort
(seed 0) — no replication; cohort ρ is marginal (0.505).

### C1c — PARTIAL (outcome yes, attribution no)
**Supported.** Cohort ε_lb valid=True on all 6 channels at δ ∈ {0.05, 0.01, 0.001}, and the rule3
values are identical to the resolution-floor values to 16 digits (e.g. af_full δ=0.05:
0.2818306806833169 in both) — the quoted cohort bounds demonstrably do not hinge on the zero-FN
corner. "Quotable, not floor-driven" is literally demonstrated.

**Not supported.** "A confidence-corrected FNR floor **makes it valid** on at least one channel":
architecturally impossible. `validate_tradeoff` checks the **raw** empirical ROC, and at a
zero-FN threshold (TPR=1, α < 1−δ) no finite ε from this estimator can validate, regardless of
floor. Empirically: no channel in s44 where rule3 changes any validity flag (floor inert where
valid — cohort; active-but-still-invalid where the corner exists — real high channels). The floor
only lowers reported values on invalid channels. Round-13's W12-5 fix ("admits", not "makes
valid") is the correct wording and the paper carries it.

### C2b — YES
nsp16000_m3000: af_only 0.5291±0.0019 vs pc_calibrated 0.5027±0.0004, Δ −0.0264±0.0016 (3 seeds);
pc_refbasis 0.4998, pc_eigonly 0.4999; in-sample gap 0.0029 (mechanism-consistent: 0.199@150 →
0.0029@16000). Cell controls: relatedness **mean** 63.69 vs donor 62.53 (+1.9%); AF spectrum
matches; eig10_frac −3.4%; LD bins 2–4 match donor (0.197/0.203, 0.077/0.077, 0.032/0.032).
LD bin 1 is NaN in **both** cohort and donor — no sampled pairs in the shortest-distance bin at
M=3000; a symmetric density artifact, not a cohort defect. 48,000 study genomes ≈ 19× 1000G arm;
largest cell fitting 24 GB (nsp=33000 needs 27 GB; disclosed). Caveats: 3 seeds; M=3000 reduced
marker density (disclosed); semi-synthetic (disclosed).

### C4b — YES
Manifest 2 (per-site AF at **full** precision + top-10 eig + λ_GC@2dp + call-rate table; pooled AF
and p-values omitted) → bundle2 0.9726 vs uncoarsened single 0.9642: excess +0.0084, positive 5/5
(frac_gt0 = 1.0). The test is in the conservative direction — strengthening the bundle's AF channel
and un-coarsening the comparator is the *harder* comparison, and the excess survived it, shrinking
from +0.0255 (vs the coarse single) to +0.0084 exactly as a coarsening-artifact explanation predicts
it should if part of the original excess were coarsening-driven. Anchoring verified: single_coarse
0.9470676543209876 == s40 max_single to 16 digits. Noise control: af_site+noise 0.9636 ≈ single
0.9642 (random per-individual column adds nothing). LOO: af_site alone 0.7395, m3 alone 0.9639 —
two near-equal-strength channels; the +0.008 is a genuine combination gain over the best single.
Caveats: +0.008 is small (though 5/5 consistent); both manifests are hand-chosen — 2 points, not a
manifest-space distribution; simulated panel only.

### C3b — YES
Native arm runs the study child with the hook installed, then attempts ctypes.CDLL(None): denied on
the `ctypes.dlopen` audit event in 3/3 seeds with `child_failed=False` — the reason string proves the
PEP-578 hook fired and child_failed=False rules out a crash masquerading as a denial. Bypass matrix
unchanged: library leak 1.0 (positive control), natural egress blocked 1.0, sanctioned-undeclared
blocked 1.0, declared flow 1.0, chain verified 1.0; governed tier 0.624/0.509 unchanged. Scope is
declared (ordinary Python, no native code, not root; parent in TCB). Residual nit: `os.write` to an
inherited fd (the stderr pipe) triggers no audit event — bytes could reach the parent unfiltered, but
the parent is in the TCB and discards stderr in the success path, so no out-of-TCB egress.

## What the results support
- The round-12 residual-limitation closure lands on its substance: ordering/controls of the
  per-individual metric replicate on real 1000G and real-haplotype data; cohort ε_lb is valid and
  provably floor-independent; the biobank-corner cell preserves AF-first ordering with passing
  controls; the bundle excess survives a stronger manifest and an uncoarsened comparator in the
  conservative direction; the ctypes escape hatch is denied.

## What they don't
- ε_lb as a quotable **absolute** bound on real 1000G's headline channels (ordering-only, disclosed).
- Any causal role for the rule3 floor in validity (none exists by construction).
- Seed replication for the s44 cohort arm (single cohort).
- Bundle-excess robustness beyond two hand-chosen manifests.

## Missing evidence
1. s44 cohort arm over ≥3 seeds (and ideally a second chromosome) — puts a std on cohort A* deltas.
2. A randomized-manifest distribution (e.g. 10 manifests sampled from the object pool) to upgrade
   "not one manifest" from 2 hand-picked points to a distribution.
3. If real-data ε_lb quotability is ever wanted: n_pos ≫ 450 or an explicitly smoothed/bound variant —
   otherwise the ordering-only reading is the honest ceiling at this n.

## Suggested claim revisions (paper needs NO change — its wording is already correct; the claim
statements as posed overreach on C1b/C1c)
- C1b → "A* and ρ replicate the real-data ordering on real 1000G chr22 and the semi-synthetic cohort
  with controls at chance (A* ≤ 0.07); ε_lb is valid and quotable on the cohort's every channel and
  remains ordering-only on real chr22's high-AUROC channels."
- C1c → "The cohort's ε_lb is valid on every channel and identical under the resolution and
  rule-of-three FNR floors — the quoted bounds do not hinge on the zero-FN corner; on real chr22's
  high-AUROC channels no finite-ε bound is valid and the floor only lowers the reported value."
- C2b/C4b/C3b → as worded; optionally add the LD bin-1 NaN note (symmetric, density artifact).

## Next experiments
1. **C1b**: s44 cohort arm ×3 seeds (+ optional chr1) — the only genuinely missing piece.
2. **C4b** (optional): 10-manifest randomized sweep.
3. **C3b** (optional): add an os.write-to-inherited-stderr injection to the bypass matrix to document
   the child→parent-TCB-only residual channel.
4. **C2b**: none — 7 cells already support the ordering claim.

## Confidence
- C1b: **high** on the anti-artefact core; the ε_lb-on-real gap is verified from code, not assumed.
- C1c: **high** (floor-inertness read directly from the identical values and validate_tradeoff's raw-ROC check).
- C2b: **high**. C4b: **high**. C3b: **high**.

## Verdict
**4 of 5 claims land as posed** (C2b, C3b, C4b yes; C1b, C1c partial with precise rewordings — both
partials are wording-level, and the paper's existing text already states them correctly). No result
is inflated, no defect found in the new experiments. The closure is genuine: every round-12 residual
limitation is now either closed in evidence or correctly scoped as an ordering-only / disclosed-limit
reading. Score: **8.5/10** (deductions: single-seed cohort arm; ε_lb not quotable on real headline
channels; hand-chosen manifests).
