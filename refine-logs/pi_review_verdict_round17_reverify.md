# Round-17 Review — independent re-verification pass (run_20261008_0d685cee)

**Reviewer**: pi. The paper is byte-identical to the round-16 state (md5 9f6cebaf…, 17 pp,
0 errors). This pass replaced eyeball-checking with a mechanical sweep: 122 quoted numbers from
the paper's abstract/Contributions/§governance/§scale and the still-present round-11/12
sections were checked against their JSON keys by script.
**Date**: 2026-10-08

## Mechanical sweep result: 120/122 pass, 2 one-character mismatches

**W17-1 (minor — number transcription).** §governance (s44 sentence): "on the cohort AF
$0.14$ leads PC $0.04$". The JSON (`s44_individual_real.json`, cohort_nsp2000): af_full
A* = **0.1328** → "0.13" at the paper's 2dp convention; PC = 0.0398 → "0.04" ✓. The results
doc says "af_full .591/.135" — the AUC .591 ✓, the A* .135 should be .133. The chain is
JSON 0.1328 → doc ".135" → paper "0.14" (double-rounded). The ordering claim (AF leads PC:
0.590 vs 0.520 AUROC; A* 0.133 vs 0.040) is unaffected. Reviewer note: this value passed my
own visual checks in rounds 12–16 (I anchored on the doc's ".135"); the mechanical sweep
caught it.
**Fix**: "0.14" → "0.13" (paper); ".135" → ".133" (results doc).

**W17-2 (nit — boundary rounding).** §governance: "the eigenvalue-only control … at
$0.55$--$0.62$". s49's six cells: eig min = **0.5447** (chr1:GBR) → "0.54" at 2dp; the quoted
lower bound 0.55 excludes the actual minimum by 0.005 (a double-rounding artifact:
0.5447 → 0.545 → 0.55). The upper bound 0.62 correctly covers s50 chr1:YRI 0.618.
**Fix**: "0.55" → "0.54" in that range's lower bound.

Everything else verified — highlights: s39 native denial 3/3 with reason `ctypes.dlopen` and
child_completed; s38 (0.61/5.0, 0.55/4.1, 0.19/1.1, controls ≤0.04/≤0.24, A* direction); s23
(0.845/0.882→0.622/0.510); s24 (0.941/0.882, 5/5); s40 (0.962/0.947/+0.015/5-of-5/LOO
0.076/noise 0.839); s46 (0.973/0.964/+0.008/5-of-5/+0.026); s26 (0.866/0.486/0.463/0.653);
s41 (planted 0.715/0.729/0.702; real 0.456/0.479/0.424 vs AF 0.866/0.864/0.838); s42 (ρ sweep
and bonus sweep exact); s37 (12-cell grid, nsp16000 48,000 genomes, relatedness +10.6% at
nsp150); s48 (12/12 single-pop gaps +0.050…+0.134 at n_study 28–36; 2-pop +0.112…+0.172;
3/5pop −0.005/−0.030 at 450/715; Fst ladder 0.004–0.087); scale-matched (+0.161@Fst 0.0997
n=60; +0.111@0.0885 n=100; controls <0.57); s49 (top-3 0.881–0.959, full 0.994–1.000, refb
0.552–0.595, oos 0.551–0.581, samehalf 0.998–1.000, chr1/2 only); s50 (within +0.051…+0.086
all positive; across-population 3/5-pop cells −0.039…+0.023 — the paper's range is correct,
my sweep's first "across" filter wrongly included the 2-pop mixes, which sit on the PC>AF
side); s21 brackets [0.027,0.052]/[0.052,0.077]/[0.103,0.151] recomputed independently —
exact; s25 (ΔE +0.24; null family ≈−0.248, paper's "≈−0.24" fair); s10 (0.879/0.832/+0.046/
+0.103/25-of-25); s16 (chr1/2 AF/PC). Compile 17 pp, 0 errors.

## 1. Score: 8.5 / 10 (unchanged — the two findings are ulp-level transcription slips)

## 2. Remaining LIMITATIONS (unchanged from round 16; blocker vs disclosed scope bound)

| # | Limitation | Class |
|---|---|---|
| L1 | No real biobank-scale real-genome cohort: the size-resolved law's large-n/low-structure quadrant is simulator-only; real within-population cohorts capped at ≈28–36 members; real large-n side is mixed 3/5-pop panels + the semi-synthetic mosaic arm (whether AF beats PC at Fst≈0 given large n is undetermined on real data) | **Disclosed scope bound** (semi-synthetic labeled; real biobank cut as unobtainable) — the largest external-validity limit, not a blocker |
| L2 | No per-split dispersion on the real round-13 gaps (5-seed means; near-crossing signs ±0.02) | Partially disclosed ("sign unstable"); cheap fix (per-split rows via one-line `_agg` change) — not a blocker |
| L3 | PC-arm orientation on the real ladder uses evaluation labels (`_cal`), inflating near-chance values ~±0.016; headline cells unaffected | Disclosed — not a blocker |
| L4 | Real-data M axis is one thinning contrast; s48/s49/s50 at one density each | Disclosed — not a blocker |
| L5 | f-DP ε_lb resolution-limited (ordering-only at zero-FN corners; flags in JSON) | Disclosed — not a blocker |
| L6 | Governance = plaintext release-object prototype, stated TCB (native boundary now tested 3/3) | Disclosed — not a blocker |
| L7 | Missingness channel constructed-by-design (ρ=0 null at 0.502) | Disclosed — not a blocker |
| L8 | 3 seeds/cell in s37/s41/s42 | Minor — not a blocker |

Plus the two transcription fixes above (W17-1, W17-2) — the only items needing edits.

## 3. Text-vs-JSON inconsistencies

Exactly two, both one-character (W17-1, W17-2). All other 120 checked numbers match their
JSON keys at the quoted precision.

## 4. Minimum fixes

- W17-1: paper "0.14" → "0.13"; results doc ".135" → ".133".
- W17-2: paper eig-range lower bound "0.55" → "0.54".
- (Recompile; nothing else.)

## Verdict

**Almost — two one-character transcription fixes away.** The paper is otherwise byte-identical
to the round-16 ready state and every conclusion survives the sweep; but a quoted number that
does not match its JSON key (0.14 vs 0.1328) is the exact defect class this loop exists to
catch, so the fixes should land before submission. After they do, the paper is submission-ready
at 8.5 with the limitation set above (all disclosed scope bounds; no blockers).
