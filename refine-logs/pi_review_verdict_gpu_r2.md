# Re-Review (GPU/scale round, Round 2) — verdict

Read `REVIEWER_MEMORY.md` (Round 1 suspicions), re-read all 9 claimed-changed files in full
(`gpu_round_summary.py`, `s29_negation.py`, `s29_paired_slow.py`, `s29_ladder.py`,
`s27_gpu_scale.py`, `s36_phased_ibd.py`, `make_figures.py`, `main.tex` sec:scale + Limitations,
`GPU_ROUND_RESULTS.md`, `PAPER_INDEX.md`), and **independently re-pooled the s29 per-seed deltas
with my own code** (not the author's aggregator). Recompiled the paper fresh (14 pp, 0 errors,
0 undefined refs). Rendered and inspected `fig8_scale.png` directly.

## Fix verification (checked against artifacts, not the change list)

| Round-1 item | Status |
|---|---|
| MAJOR: 600-seed sweep was one 200-seed shard (dict.update) | **FIXED & VERIFIED** — `merge_paired` concatenates `per_seed_delta`; my independent pooling gives +0.0943 [+0.0925,+0.0962] P=1.0 (n=600), +0.0088 [+0.0083,+0.0094] P=0.943 (n=600), −0.0136 P=0.210 (n=200) — all match the paper. |
| MAJOR−: "saturates near k≈8–10" unsupported | **FIXED** — increments [0.132,0.066,0.055,0.043,0.037,0.028,0.019,0.018,0.014] match the quoted "+.13,+.07,+.06,+.04,+.04,+.03,+.02,+.02,+.01"; strictly decaying; "still climbing slightly at k=10" is honest (+0.014). |
| eig-only range wrong (0.508–0.533) | **FIXED** — real chr22/1/2 = 0.503/0.533/0.521, msprime 0.508, sim 0.510 → "0.50–0.53" correct; biobank grid "0.497–0.511" matches 0.4965–0.5115. |
| "loadings at every scale" overreach | **FIXED** — `pc_eigonly` arm now exists in all 12 biobank cells (real control, orientation-corrected, GPU-scored); text qualified to "every tested corpus and scale". |
| Undisclosed gauss approximation | **FIXED & ADEQUATE** — Limitations names the Gaussian diffusion approximation, scopes it to biobank-scale cells, states O(Fst²), and I verified `common.py` (used by the 600-seed sweep) uses the exact `rng.beta` BN draw. Does NOT materially weaken the headline: the paired-negation cells are exact-Beta, the boundary is a large-effect phenomenon, and msprime + real-1000G agree externally. |
| Limitations unbenchmarked GPU clause + dangling ref | **FIXED** — now "we did not benchmark a GPU parser"; parser numbers verify (pigz pipe 2471.6/s → "2.5k", gzip 3290/s → "3.3k"). |
| Consolidation defects (S36 empty, S28 "1e4", S29 "600/cell") | **FIXED** — S36 section populated; S28 header "20000"; S29 table shows true n_seeds (600 fast / 200 slow). |
| Fig 8 caption mismatches | **(c) FIXED** (no "log scale"); **(a) NOT FIXED — see W1**. |
| s31 per-pop disclosure | **Sidestepped, acceptable** — quoted numbers are the 5-pop arms (per=100, so "matched at 100" is true of the quoted comparison); the chr22 3-pop arm (per=150) is no longer quoted (nit). |
| s36 non-member −1 dummy | **FIXED** — `self_blocks=None` scores non-members against the full released set; members lose both own haps; symmetric and correct. Numbers verify (real 0.486, planted 0.671). |

Also verified: s30 density numbers (20×, PC flat 0.66–0.67, AF rises); s35 frontier (+0.008 at
M=5k, −0.02 at M=10⁴, −0.03 at M=5·10⁴, −0.03→−0.17 widening for nsp2000); ref-basis collapse
0.4854–0.5041 in all 12 cells; abstract "eigenvalue-only release is at chance" now backed;
`PAPER_INDEX.md` artifact map accurate; compile clean.

## 1. Score: **8/10** for a top venue

Both round-1 MAJORs are genuinely fixed and independently verified; the eig-only control is a real
new measurement backing the governance-critical loadings claim. Every headline number in sec:scale
reproduces from pooled raw data. It stays below 8.5 on a recurring caption/provenance slip pattern —
not on any wrong measurement.

## 2. Remaining critical weaknesses (ranked)

1. **W1 (MAJOR−) — Fig 8(a) caption is false about its own panel.** Caption: "the
   reference-fitted-basis control gap (dashed) sits at zero everywhere—i.e. that control channel is
   at chance." The plotted dashed series is `pc_refbasis − af_only`, which spans **−0.25 to −0.51**
   (confirmed in `make_figures.py` and in the rendered PNG — the dashed curves sit far below the
   zero line). The control channel *is* at chance in AUROC (0.485–0.504), but the plotted series is
   not zero; the caption contradicts the figure a reader is looking at. Round-1 item 8 was claimed
   fixed but only the series *label* changed.
2. **W2 — "pooling all six shards of a paired sweep" is wrong.** The artifact trail holds **3**
   `s29_negation_shard*.json` files (200 seeds each) plus 1 slow run. Same provenance-mislabel
   pattern as round 1.
3. **W3 — seed count for the n=1600 cell is unstated and smaller.** The sentence's "(600 seeds)"
   attaches to +0.094; the M=8000 cell does have 600, but the n=1600 (and Fst=0.24) cells have
   **200** seeds each (from `s29_negation_slow.json`). A careful reader assumes 600 throughout.
4. **W4 — Limitations says "+0.095"**, sec:scale says "+0.094"; actual is +0.0943. Round-1 fix
   updated sec:scale but missed Limitations.
5. **W5 (MINOR) — in-sample-gap endpoints cherry-pick columns.** "from +0.51 at n_sp=500 to +0.08
   at n_sp=8000" mixes 0.513 (f0.04/M100k) with 0.082 (f0.08/M50k). Per-column the shrinkage is
   monotone (verified), but the n_sp=8000 gaps span 0.082–0.148, so "+0.08" quotes the floor.
6. **W6 (MINOR, latent) — shard-winner inconsistency moved to s27.** `s27` shard2 re-ran 2 cells
   already in shard0 (different RNG draws; `pc_eigonly` absent in shard2). `gpu_round_summary`
   takes first-wins, `make_figures.merge_shards` takes last-wins → figure and summary disagree by
   ~±0.004 on those cells. No conclusion changes, but it is the same dedup-policy bug class as
   round 1's s29 defect, now in the figure path.
7. **NITs** — (a) biobank grid runs 3 seeds/cell (plan said 8) and the paper never states the grid
   seed count; (b) `s29_negation.py` docstring says the eig-only control uses "the same spectrum"
   but the code draws a plain random orthonormal basis (mislabeled comment, not a wrong control);
   (c) s31 chr22 3-pop arm remains per=150 (unquoted).

## 3. Minimum fixes

- **W1**: one-line caption edit — e.g. "the dashed series shows ref-basis−AF; the control channel
  itself sits at chance (0.485–0.504)" — or plot `refbasis−0.5` with the zero line. No re-render
  needed if the caption route is chosen.
- **W2/W3**: reword to "pooling three 200-seed shards (fast cells) plus a 200-seed run of the two
  slow cells (600 seeds per fast cell, 200 for n=1600 and Fst=0.24)".
- **W4**: "+0.095" → "+0.094" in Limitations.
- **W5**: "from +0.49–0.51 to +0.08–0.15" or quote one matched column.
- **W6**: `setdefault` in `make_figures.merge_shards`, or delete the 2 duplicate cells from shard2
  (then regenerate fig8).
- **NITs**: one sentence stating the grid seed count; fix the docstring wording.

All are text/one-line-code fixes; no re-experimentation required.

## 4. READY? **Almost**

Submission-plausible after one ~30-minute text/plot-fix pass. The measurement core is sound and now
honestly aggregated; the Fig 8(a) caption contradiction is the only defect a venue reviewer would
catch at a glance.

## Memory update

- **Round-1 suspicions resolved?** Pooling (dict.update): fixed and independently verified.
  k-ladder saturation: honestly reworded. Biobank eig-only: real arm, all 12 cells. Gauss: verified
  adequate (small-scale cells confirmed exact-Beta in common.py).
- **New suspicions**: (1) fix-lists over-claim — "Fig 8(a) now describes the panel honestly" was
  false; verify captions against the RENDERED panel, not the intent; (2) shard-winner dedup
  inconsistency persists in the figure path (first-wins vs last-wins on s27 shard2 duplicates);
  (3) prose seed/shard counts still drift from the artifact trail ("six shards", implicit 600 for
  200-seed cells).
- **Patterns**: caption/figure mismatch class survived a claimed fix; always render and inspect the
  figure, not just the plotting code. Provenance counts in prose need a per-number audit against
  file inventories.
