# Cross-model review bundle — round 3 (paper-level review)

You are a PoPETs/USENIX reviewer. Review the accompanying DRAFT PAPER
(`paper/main.tex`, 4 pages) and the evidence it cites. The executor implemented
your round-2 must-dos. Judge whether the paper is now submission-plausible and
name the single biggest remaining weakness. Do NOT propose adding modules;
rank and judge what exists.

## What changed since round 2
1. **S3 implemented non-degenerately** (`code/s3b_paradox.py`): the analyst's
   stopping rule is COMPUTED from the revealed spectrum (two families: scree
   threshold, largest-gap), quality is MEASURED by re-running the association
   pipeline at the chosen PC count (calibration + recovery), probes are
   emergent. Sanitized arm pre-provides a coarse diagnostic block.
   Result: probes 3→0 (A) and 3→0 (B), equal/better measured quality.
2. **S4 re-run with the disjoint split** (`s4_msprime.py`): full-AF AUROC 0.64,
   1dp 0.65 (was 0.72/0.63 with the biased neg set).
3. **S1c**: frontier restored to 5 seeds + per-ancestor breakdown, disjoint
   negatives (panels reused across tiers for speed).
4. **S6**: dimensionality sweep for the PC-leak magnitude.
5. **S2c**: mechanism formalised (reviewer's own diagnostic promoted to code).

## Current headline numbers (all in results/raw/, verified)
- Frontier (5 seeds, disjoint neg): F2 0.582 / F4(1dp) 0.527 / F3(2dp) 0.714 /
  F5 0.720 / F5-per-site 0.858; ancestry control 0.494; per-ancestor uniform.
- PC-spectrum leak (5 seeds): base 0.753 → +PC 0.861 (Δ+0.108); pcdist-only 0.79
  (> base); pcalign-only ≈chance; ref-fit control Δ≈0 (T2); study-mean centering
  Δ intact (T3); eig ratio ≈1.00.
- Site×precision boundary; dimension sweep Δ from +0.03 (m/n=60) to +0.12.
- msprime cross-check 0.64/0.65.
- S3b: probes 3→0 both families, quality parity.

## Draft paper claims
- Headline: "the diagnostics (PC-spectrum) release is the dominant leak" — PC
  projection attack (0.79) beats full AF release (0.75); combined 0.86.
- Self-reference-bias correction (≈0.2 AUROC) as a methodological contribution.
- Site×precision boundary as actionable.
- S3 framed as a model-level conjecture only.

## Your task
1. Is the paper's headline claim (_diagnostics release dominates_) supported by
   the evidence? SUPPORTED / PARTIAL / NOT.
2. Is the self-reference-bias correction a real contribution or a
   self-inflicted-repair?
3. Single biggest remaining weakness for a PoPETs submission.
4. What is the reviewer's most likely rejection reason, and is it fatal or
   fixable within the current evidence?
5. Verdict: submission-plausible as a measurement paper, or not yet?
