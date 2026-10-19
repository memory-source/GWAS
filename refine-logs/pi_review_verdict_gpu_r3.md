# Re-Review (GPU/scale round, Round 3) — verdict

Read all 5 claimed-changed files, re-verified every number against raw JSON, **rendered and
inspected `fig8_scale.png` directly**, re-pooled nothing new (s29 shards unchanged, timestamps
05:19–05:41), recompiled fresh (14 pp, 0 errors, 0 undefined refs).

## Fix verification (checked against artifacts, not the change list)

| Round-2 item | Status |
|---|---|
| W1: Fig 8(a) caption false | **FIXED & VERIFIED** — plot now draws `pc_refbasis − 0.5` (code confirmed); caption "Dashed: the reference-fitted-basis PC AUROC minus 0.5 (distance from chance), which sits on the zero line in every cell" is TRUE: refbasis ∈ [0.4854, 0.5041] → refbasis−0.5 ∈ [−0.015, +0.004]; the rendered PNG shows the dashed curves hugging the zero line. Caption and panel now agree. |
| W2: "six shards" | **FIXED** — now "three 200-seed shards + one 200-seed run for the slow cell, 600 seeds on each fast cell"; matches the artifact trail (3 `s29_negation_shard*.json` + `s29_negation_slow.json`; 600 seeds per fast cell, 200 per slow cell). |
| W3: n=1600 seed count | **FIXED** — "(200 seeds)" now explicit next to the −0.014/P=0.21 number. |
| W4: Limitations "+0.095" | **FIXED** — no "0.095" anywhere in main.tex; Limitations and sec:scale both say "+0.094" (= actual +0.0943). |
| W5: gap endpoints mixed columns | **IMPROVED** — now per-Fst columns. Fst=0.04 "+0.51→+0.15" is correct (0.513→0.148, matched M100k). Fst=0.08 "+0.50→+0.08" quotes the M50k cell (0.082) while the matched-M100k endpoint is +0.11 (0.110) — minor residual (nit below). |
| W6: shard-winner flip | **RESOLVED IN EFFECT** — shard2 re-run now contains only `nsp2000_m100000_f{0.04,0.08}` with values **byte-identical to shard0** (all 5 fields SAME, verified). First-wins vs last-wins is moot. NOTE: the claimed mechanism ("make_figures uses setdefault, first-wins") is **not in the code** — `merge_shards` still uses `out.setdefault(key, {}).update(d[key])`, i.e. last-wins at the cell level; the no-flip outcome holds only because the overlapping values are identical. Latent policy inconsistency remains (no current impact). |
| NIT: s29 docstring | **FIXED** — "random orthonormal basis, no fitted directions". |
| NIT: GPU_ROUND_RESULTS / summary | **Regenerated correctly** — S27 all 12 cells with eig-only (0.497–0.511), S29 pooled table with true seed counts, S36 populated. |

Also re-verified: eig-only range 0.4965–0.5115 over all 12 cells (no missing); ladder increments
unchanged and honestly described; per-Fst in-sample-gap monotonicity true in every column;
"+0.094 [+0.092,+0.096] P=1.0 (600 seeds)" still matches my independent pooling; compile clean.

## 1. Score: **8.5/10** for a top venue

Three rounds of adversarial review with independent recomputation have converged: 2 MAJORs (R1) →
1 MAJOR− (R2) → 0 (R3). Every headline number traces to pooled raw data; captions now match the
rendered panels; provenance counts match the artifact trail; the gauss approximation is disclosed
and verified immaterial to the headline; the loadings claim is backed at every tested scale. It
stays below 9 on residual polish items only.

## 2. Remaining weaknesses (all nit-level, ranked)

1. **Fst=0.08 gap endpoint "+0.08" uses a different M than the Fst=0.04 endpoint.** Fst=0.04
   quotes the matched-M100k endpoint (+0.15 = 0.148); Fst=0.08 quotes the M50k cell (0.082) while
   the matched-M100k endpoint is +0.11 (0.110). A recomputing reviewer would notice the
   inconsistency. No conclusion change (shrink from ~+0.5 to ~+0.1 either way; per-column
   monotonicity verified).
2. **Latent dedup-policy inconsistency**: `gpu_round_summary.merge_s27` first-wins vs
   `make_figures.merge_shards` last-wins. Harmless now (overlapping values byte-identical) but the
   round-3 claim that make_figures "uses setdefault (first-wins)" is not what the code does.
3. **Claim-list drift (not a paper defect)**: the change report describes shard2 as containing
   "nsp500_m1e6, nsp2000_m100k f0.04/f0.08"; actual shard2 holds only `nsp2000_m100k` × 2. The
   paper never describes shard2's contents, so the artifact trail and paper are unaffected.
4. **Carried nits**: "one 200-seed run for the slow cell" (singular — the slow run covers two
   cells, n1600 quoted and f0.24 unquoted); biobank grid seed count (3 seeds/cell) still unstated
   in the paper; s31 chr22 3-pop arm per=150 (unquoted).

## 3. Minimum fixes (all optional polish)

1. "+0.08" → "+0.11" (or "+0.08–0.11") in the Fst=0.08 clause — one word.
2. Unify `merge_shards` to `setdefault` for hygiene — one line, no output change.
3. Optional: one clause stating the grid runs 3 seeds/cell; "the slow cell" → "the slow cells".

## 4. READY? **Yes**

Submission-plausible now. Every measured number verifies against pooled raw data; captions match
rendered panels; provenance counts match the artifact trail; disclosures are adequate; the
controls are real at every tested scale. The remaining items are copy-edit level and would not
block acceptance; the two one-word/one-line polish items are recommended before submission but not
required for correctness.

## Memory update

- **Round-2 suspicions resolved?** All: W1 caption verified against the rendered PNG (dashed now
  refbasis−0.5, sits at ~0); W2–W4 fixed; W6 resolved via value-identical shards.
- **New suspicions**: (1) fix-claim mechanisms can be inaccurate even when the outcome is real —
  the W6 "setdefault" claim was not in the code; verify the mechanism, not just the outcome;
  (2) change-report cell inventories can drift (shard2 contents misdescribed) — check file
  contents, not the report.
- **Pattern confirmed**: the loop converged 2 MAJORs → 1 MAJOR− → 0 across three rounds with
  independent recomputation at each step; the recurring failure class (caption/figure mismatch,
  provenance mislabeling) was eliminated only after captions were checked against rendered output.
