# Senior Review — GPU/scale round, Round 1/4 (run_20261006_350bbd97)

Reviewed directly from artifacts: `paper/main.tex` (13 pp) in full; all 10 listed
scripts plus `s32_phased_ibd.py`, `s36_phased_selfexcl.py`, `gpu_round_summary.py`,
`make_figures.py`, `common_gpu.py`; all 13 raw JSONs verbatim; recompiled fresh
(pdflatex ×2). Script-matched every number quoted in Section `sec:scale` against
`results/raw/*.json`, and re-pooled the S29 shards independently (the
consolidated summary provably does not pool them — see W1).

## 1. Score: 7.5/10

Content-wise this is the strongest round yet: the boundary grid is a real
experiment with clean controls, the stratified permutation null is correctly
implemented, and every quoted grid/frontier/density/5-pop/IBD number traces to
raw JSON. It is held below "ready" by one systematic provenance mislabeling
(W1), two overstatements the data do not support (W2, W4), one undisclosed
simulator approximation (W5), and one unsupported Limitations clause with a
dangling cross-reference (W6).

## What holds up (verified, not trusted)

- **Boundary grid (S27/S35)** — every quoted number matches raw: +0.008 at
  nsp500/m5000; −0.0221/−0.0304 at m1e4/m5e4; nsp2000 −0.0329→−0.1736 monotone
  in M; nsp8000 AF 0.751–0.837 vs PC 0.585–0.646; ref-basis control 0.4854–0.5041
  in **all 12 cells**; in-sample gap 0.513→0.082 monotone across study size.
  Disjoint splits, train-half orientation calibration, out-of-fold logistic all
  verified in code. Free replication: shard2 re-ran 2 of shard0's cells and
  reproduces within ±0.008 (delta +0.0041 max).
- **Permutation nulls (S28)** — NPERM=20000 (matches the paper's "2·10⁴");
  stratification permutes within ancestry (correct null); p = 4.9998e-5 < 5e-5 in
  all 5 cells, so "p<5·10⁻⁵" is literally true.
- **Paired sweep (S29)** — genuinely 600 seeds/cell by design (3 shards × 200),
  paired per-seed deltas, P(Δ>0)=1.0 holds pooled (600/600). Point estimates
  survive pooling: +0.0943/+0.0088/−0.0138 (see W1 for the CI/P provenance).
- **Ladder + eig-only control** — monotone in k on all five corpora; eig-only
  0.503–0.533 at chance everywhere tested (see W3 for the quoted range).
- **Density + 5-pop (S30/S31)** — PC flat 0.656–0.667 while AF rises at every
  tier; gap widens on all 3 chromosomes; 5-pop preserves the ordering in all 4
  panels (Δ −0.033…−0.186, all negative).
- **Phased IBD (s36)** — the correction is real: members' own two haplotypes are
  removed from the released set (the s32 self-match flaw is honestly documented
  in the s36 docstring and visible in s32's raw JSON error + 0.599-at-all-
  thresholds output). Raw: real 0.486 ± 0.033 (5 splits 0.441–0.517), planted
  control 0.671 ± 0.025 — matches the paper's 0.486/0.671.
- **Negative results honestly recorded in artifacts**: s32's failed positive
  control; pigz *slower* than the Python parse (2471 vs 3290 var/s); s33's
  MLP genuinely computed (independent stds), though never cited by the paper.
- **Harness**: TF32 off, host-resident int8, exact Gram `eigh`, streamed blocks —
  as claimed. Compile: 13 pp, 0 errors, 0 undefined refs, 4 cosmetic overfull.

## 2. Remaining critical weaknesses (ranked)

**W1 (major). "600-seed paired sweep" numbers are actually one 200-seed shard.**
`gpu_round_summary.py` merges the three s29 shards with `dict.update()`, which
overwrites per-cell — last shard wins — so the consolidated S29 table (labeled
"600 seeds/cell", n column = 200) holds only shard2 (seeds 400–599). The paper
quotes exactly those shard2 numbers as 600-seed in two places: sec:scale
("+0.095 [+0.092,+0.098], P=1.0 … +0.009 [+0.008,+0.010], P=0.95") and the
Limitations ("decays from +0.095 to +0.009"). Pooled over all 600: +0.0943, CI
≈ [+0.092,+0.096] (the quoted hi 0.098 is *wider* than the 600-seed CI), and
P(Δ>0)=0.943 at M=8000 (quoted 0.95). Conclusions identical; provenance wrong in
4 quoted numbers. Trust-relevant: a reviewer who pools the shards (as I did)
finds the paper's CI is not the claimed experiment's CI.
*Minimum fix*: pool the shards in the summary script (store per-seed deltas, not
per-shard aggregates), regenerate `GPU_ROUND_RESULTS.md`, update the four quoted
numbers. No GPU rerun needed — the shards contain all 600 seeds' aggregates;
ideally also emit per-seed deltas for an exact pooled CI.

**W2 (major−). "Saturates near k≈8–10" is not supported.** Increments at
k=9→10 are +0.014 (sim), +0.015 (chr22), +0.016 (chr1/chr2), +0.011 (msprime) —
as large as k=8→9 in 3 of 5 corpora. The ladder is still climbing at k=10.
Same word appears in the Fig 8 caption ("saturating").
*Minimum fix*: soften to "increments decay monotonically" (text + caption), or
extend one ladder to k=15–20 to demonstrate an actual plateau.

**W3 (minor). Eig-only range slip.** Paper: "(0.508–0.533)" for
sim/real-chr22/chr1/chr2/msprime. chr22's eig-only is **0.503** — below the
stated lower bound. True range: 0.503–0.533.
*Minimum fix*: one number ("0.50–0.53" or "0.503–0.533").

**W4 (minor+). "The leak is in the loadings, at every scale" — overreach.**
The k-ladder and eig-only control run at small-n corpora only (sim 360 study,
real 450–750, msprime); there is **no eig-only arm at biobank scale** (S27 has
the ref-basis control, which is a different control — fitted loadings, not
eigenvalue-only). The sentence sits in the section whose whole point is scale.
*Minimum fix*: reword to "at every tested corpus and operating point", or add
one eig-only arm at a single biobank cell (cheap: reuse `gram_svd_host`'s
eigenvalues with a random orthonormal basis, as s29 does).

**W5 (minor+). Undisclosed simulator approximation + stale seed inventory.**
S27/S35 call `simulate_population_host(..., pct="gauss")` — a Gaussian
diffusion approximation of the Balding–Nichols draw (exact to O(Fst²) per the
docstring), while §controlled states "Balding–Nichols" and sec:scale discloses
nothing. Also the §controlled seed inventory covers only the *old* experiments;
the new ones are absent (s27: 3 seeds/cell, s28: 3, s29: 600, s30/s31: 6,
s35: 6, s36: 5 splits). The plan's "8 seeds/cell" for B1 and "1e5 permutations"
for B3 were silently not met (paper says 2·10⁴, which *does* match code).
*Minimum fix*: one disclosure sentence in sec:scale ("the biobank-scale grid
uses a Gaussian diffusion approximation of the BN allele-frequency draw, exact
to O(Fst²), for tractability") + extend the seed inventory by one clause.

**W6 (minor+). Unsupported Limitations clause with a dangling cross-reference.**
"a binary VCF parse is CPU-bound (a GPU decompressor is not a win;
Section~\ref{sec:scale})" — (a) sec:scale contains no parser discussion at all;
(b) no GPU path was ever benchmarked: s34 measures only python-gzip, pigz-pipe,
and pigz-decompress (all CPU). The plan's "GPU string-parsing straw man" was
never run; "not a win" is an inference, not a measurement. What *was* measured:
pigz-p16 pipe is slower than the Python parse on the first 20k variants, and
full-chr22 decompression at 16 threads takes 173.7 s vs ≈335 s projected parse.
*Minimum fix*: move the B6 sentence into sec:scale and restate only what was
measured ("decompression at 16 threads is not the bottleneck; the parse is"),
dropping the GPU claim or marking it as not benchmarked.

**W7 (minor). Consolidation defects (same root as W1).** GPU_ROUND_RESULTS.md:
S36 section is empty (the handler expects s32's schema; s36's `agg` is
top-level), so the paper's 0.486/0.671 appear nowhere in the "verbatim"
consolidated record; S28 header says "1e4" (2e4 actual); S29 table labeled
"600 seeds/cell" but holds 200; S27 "first shard wins" silently ignores
shard2's two re-runs (they reproduce, so no data harm).
*Minimum fix*: same summary-script pass as W1 + s36 schema case.

**W8 (minor). Fig 8 caption/panel mismatches.** Caption (c) says "(log scale)":
panel (c) has linear axes. Caption (a) says "the reference-basis control is at
chance everywhere": panel (a)'s dashed series plots `pc_refbasis − af_only`
(a negative gap), so the 0.5-chance level is not visible anywhere in the panel.
*Minimum fix*: add `set_yscale("log")` (optional) or fix the caption; plot the
ref-basis AUROC itself (or annotate 0.5) in panel (a).

**W9 (minor). s31's "matched" 3-pop vs 5-pop comparison is not matched for
chr22** (3-pop at per=150, 5-pop at per=100; the code comment says "matched ...
per=100" but chr22_3pop ran at 150). The paper only claims non-reversal, which
holds; the asymmetry weakens the comparison. *Minimum fix*: state the per-pop
sizes, or rerun chr22-3pop at per=100.

**W10 (nit). s36 non-member self-exclusion is a dummy.** `self_non = −1` makes
`keep[-1]=False`, removing only the *last* released haplotype for every
non-member (numpy negative indexing), while members lose exactly two. This
slightly *lowers* the real-arm AUROC (conservative for the "at chance" claim),
but the paper's "with the candidate's own two haplotypes removed" is inaccurate
for non-members. *Minimum fix*: give non-members the full released set (they
have no self rows) and note it; effect is negligible (one hap out of ~800).

**W11 (process). Undisclosed plan deviations** (paper consistent with code, so
not result defects): B1 8 seeds → 3; B3 1e5 perms → 2e4; B5 chr1-3 with
thresholds L∈{20,50,100} → chr22 only, no thresholds, 5 splits; B6's GPU straw
man never run; the plan's "Table 1 CIs" target was not delivered (Table 1
unchanged — acceptable, prior CIs verified). *Minimum fix*: a deviations
paragraph in GPU_ROUND_RESULTS.md.

## 3. Minimum fixes (consolidated, in priority order)

1. Fix the s29 merge to pool all 600 seeds; regenerate summary; update the four
   quoted numbers (W1/W7). *Text + summary script; no GPU.*
2. Reword "saturates near k≈8–10" / "saturating" (W2). *Two strings.*
3. Eig-only range → "0.50–0.53" (W3). *One number.*
4. Reword "at every scale" → "at every tested corpus/operating point" (W4), or
   add one biobank-scale eig-only arm. *One string or one cheap run.*
5. Disclose the `gauss` approximation + extend the seed inventory (W5). *Two
   sentences.*
6. Move the parser sentence into sec:scale, restate only measured facts, drop
   the unbenchmarked GPU clause (W6). *Two sentences.*
7. Fig 8: caption "(log scale)" and the ref-basis series (W8). *Caption + one
   plot line.*
8. (Advised) s31 per-pop disclosure (W9); s36 non-member handling note (W10);
   plan-deviations paragraph (W11). *Sentences.*

## 4. READY? Almost

No result is fabricated or invalidated; every checked sec:scale number matches
raw JSON, and the S29 conclusions survive independent pooling. But a submission
carrying "600-seed" numbers that are provably one shard, a saturation claim the
ladder does not show, an undisclosed simulator approximation under the new
section's headline grid, and an unbenchmarked GPU claim is not submission-ready.
Fixes 1–7 above are hours of work (one summary-script fix + ~8 text edits, one
optional cheap run); after them this is a credible 8.5 and ready.
