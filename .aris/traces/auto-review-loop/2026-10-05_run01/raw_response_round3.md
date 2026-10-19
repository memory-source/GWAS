# Senior Review — Round 4 fresh, Round 3 (confirmation after final fixes)

Reviewer: independent artifact-level audit. Verified the two remaining line edits in the
current `main.tex`, ran my own full ordering-claim sweep (not just the four lines the
change list named), confirmed the edits are text-only (no numbers touched: 0.729±0.015,
0.964, 0.990, +0.055 all intact), and recompiled (10 pp, 0 errors, 0 undefined refs, 3
cosmetic overfull; tex+pdf regenerated together 00:33).

## Fix verification

- **Discussion (line 497) — FIXED.** "…per-site frequency release is the largest
  marginal-frequency single channel. Two other release families can be larger still---a
  haplotype/LD table and a membership-correlated missingness table---which only widens the
  set of released objects a manifest must cover." Accurate against the artifacts (hap-LD
  0.990 > per-site AF 0.964; M3 up to 1.00; "can be larger still" hedged correctly).
- **§6 taxonomy intro (line 231) — FIXED.** "reported separately as the largest
  marginal-frequency single channel."
- **W5 subset naming (line 384) — FIXED (optional item applied).** "s20_extra.py, the
  latter on a 2000-marker subset" — the s15/s20 structure reconciliation is now airtight.

## Full ordering-claim sweep (all instances, my own grep)

- Abstract line 36: scoped by the canonical-set sentence ✓
- Abstract lines 42–43: "largest marginal-frequency channel (0.96); two channels outside
  this canonical set are larger still" ✓
- Contribution 1 heading (line 95): "dominant intermediate leak" — category-level thesis
  (diagnostics-as-a-class vs the final statistic), body scoped; matches the paper title ✓
- Contribution 1 body (lines 100–102): scoped ✓
- §6 intro (line 231): fixed ✓
- Line 392: "sub-dominant on real data" — correct magnitude statement ✓
- Line 412: "we therefore scope the claim: the diagnostics release is the largest
  intermediate channel on the simulated panels; on real data the PC diagnostic is…" —
  same class as the fixed claims; the category subject makes the thesis reading available
  (same reading as the contribution-1 heading), so it is defensible, but a one-word polish
  ("canonical") would remove all ambiguity. Truly optional, non-blocking.
- Discussion (lines 497–498): fixed ✓

## Regression check
All previously fixed items hold (disjoint-negative control, loadings-as-released-object,
s7 multi-split, de-confounded s18 mechanism 0.729/0.510/0.519, F/M taxonomy, seed
inventory, missingness caveat, structure reconciliation, hap-LD operating point +
family-bound clause, bibliography). Edits were text-only; no measured number changed.
Compile clean.

## 1. Score: 9/10

## 2. Remaining weaknesses (none blocking)
1. Line 412's scope sentence — one-word optional polish ("largest *canonical* intermediate
   channel" or "largest intermediate channel among the canonical diagnostics").
2. Disclosed limitations, not defects: no MPC enforcement prototype, no phased-haplotype/IBD
   attack, magnitudes model-specific — all honestly stated in Limitations.

## 3. Minimum fixes
None blocking. Optional: the line-412 one-word polish.

## 4. READY? **Yes** — submission-plausible as a measurement paper.
