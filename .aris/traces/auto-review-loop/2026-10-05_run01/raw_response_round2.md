# Senior Review — Round 4 fresh, Round 2 (confirmation after fixes)

Reviewer: independent artifact-level audit. Re-read `main.tex` (current), the rerun
`s18_msprime_pc.py` + `s18_msprime_pc.json`, `s17_scale.py` (docstring), `s19_haplotype.py`/
`s20_extra.py` (unchanged, timestamps), `make_figures.py` fig6, and recompiled (10 pp, 0
errors, 0 undefined refs, 3 cosmetic overfull).

## Fix verification (against artifacts, not the change list)

- **W2 (s18 swapped null) — FIXED, substantively.** The mechanism now uses a
  within-population ancestry-matched split (code comment documents the S9d/S14 rationale:
  splitting by msprime population would confound the null with ancestry separation).
  Per-seed swapped values are tight (0.504/0.510/0.504/0.523/0.538 → 0.5185 ± 0.0119);
  paper updated exactly: in-sample 0.729 ± 0.015, out-of-sample 0.510 ± 0.010, swapped
  0.519 ± 0.012. The structure-contaminated 0.967/0.547/0.566 numbers are gone. The
  in-sample drop 0.967 → 0.729 is the honest de-confounded number, and the mechanism
  signature (0.729 vs ≈0.51) holds. Other s18 numbers unchanged (deterministic seeds), so
  PC 0.691/AF 0.627/granularity still match. Fig. 6(b) plots no mechanism bars → not stale.
- **W1 (ordering-claim contradiction) — PARTIAL.** Abstract FIXED: canonical-set scoping
  ("among the canonical GWAS diagnostic releases (λGC, the PC spectrum, missingness and the
  condition bit)") plus the explicit sentence "two channels outside this canonical set are
  larger still—a haplotype/LD channel (0.99) and a membership-correlated missingness table
  (up to 1.00)". Contribution 1 body FIXED the same way. **BUT the Discussion (line 496)
  still reads "per-site frequency release is the largest single channel"** — unscoped, and
  contradicted two paragraphs later by the paper's own hap-LD 0.99 / M3 1.00 — and §6's
  taxonomy intro (line 231) has the same unscoped phrase ("reported separately as the
  largest single channel"). The change list claimed the Discussion said "largest
  marginal-frequency single channel"; it does not. Two one-line edits remain.
- **W3 — FIXED.** "(+0.055–+0.056)" matches s17 (0.0565 / 0.0550).
- **W4 — FIXED.** s17 docstring now matches what ran ({1200,3000}, M=10000, ld100 + one
  zero-LD cell, 8 seeds/cell); §4 seed inventory extended (multi-chrom 20 splits/chromosome,
  site-size 8 seeds/cell, msprime/haplotype 5 seeds).
- **W5 — FIXED (minimum met).** "≈0.10 on the full chromosome vs ≈0.05 simulated; chr1 and
  chr2 are higher still, ≈0.14 and ≈0.11; … computed in s15_controls.py and s20_extra.py" —
  all quoted numbers artifact-backed (s15 0.0994; s20 0.1418/0.1139); the "full chromosome"
  qualifier distinguishes s15's estimate from s20's 2000-marker subset, and s20's chr22
  value (0.0725) is not quoted, so no visible inconsistency remains. Optional polish: name
  the subset explicitly.
- **W6 — FIXED.** Operating point stated (900/900/900, 4000 markers, 5 seeds — matches
  s19's code) + release-realism clause ("family bound rather than a standard release").

## Regression check
All previously fixed items remain fixed: disjoint-negative control, loadings as the
released object, s7 multi-split numbers, F/M taxonomy, missingness caveat, labels/refs,
bibliography. Compile clean (10 pp, 0 errors, 0 undefined refs). s19/s20 untouched; their
numbers still match the paper.

## 1. Score: 8.5/10 (up from 8.0 — real improvement: de-confounded mechanism + headline scoping)

## 2. Remaining weaknesses (ranked)
1. Two residual unscoped "largest single channel" claims — Discussion (line 496, the
   material one) and §6 intro (line 231, context-mitigated). The change list over-claimed
   the Discussion fix (the recurring fix-list-vs-diff pattern; this round 5/6 items were
   delivered exactly as claimed).
2. Minor/optional: W5 reconciliation is implicit — naming the s20 2000-marker subset would
   make it airtight.

## 3. Minimum fixes
1. Line 496: "…per-site frequency release is the largest marginal-frequency single
   channel" (or "the largest tier of the release-precision frontier"). Line 231: same
   one-clause scoping. Recompile.
2. Optional: add "(2000-marker subset)" to the s20-based structure clause.

## 4. READY? Almost — one text-fix pass (two sentences) from Yes.
