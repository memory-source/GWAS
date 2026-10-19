# Round-3 Fresh Review (autonomous review loop, round 1/3 of this run) — 2026-10-04

Reviewer: independent artifact-level review. Read main.tex, all code (common.py,
real_panel.py, s9–s14 + prior scripts), all raw JSONs verbatim, CLAIMS_FROM_RESULTS.md,
review-stage/AUTO_REVIEW.md. Re-compiled the paper (8 pp, 0 errors, 0 undefined refs).
Re-ran S9c granularity for 3 seeds: **bit-identical** to stored JSON (diff 0.00e+00).

## What holds up (verified, not trusted)

- **~65 paper values script-matched against results/raw/*.json — all match within
  rounding.** Tables 1/2/3; per-ancestor 0.860/0.849/0.851; S10 paired stats (PC
  0.879±0.009 vs AF 0.832±0.021; +0.046, CI [0.037,0.055]; t/Wilcoxon/sign <1e-4;
  25/25 sign; combined +0.103, CI [0.096,0.110]); T2 +0.003; T3 +0.102 (per-seed
  0.092–0.116 → "+0.09 to +0.12"); eig ratio 0.9998; perm 0.486; |PC1| 0.501;
  S9d 0.881/0.880 own vs 0.510/0.505 other vs 0.505/0.510 swapped; S9c ladder
  0.527/0.607/0.655/0.740/0.882 (eig 0.510); S9b 0.502 vs 0.634; s9 real 0.652±0.014
  vs 0.688±0.026 (20/20 farther); s13 ladder 0.509/0.526/0.545/0.579/0.652 (eig 0.515);
  s14 0.685±0.015 vs 0.516±0.013 vs 0.516±0.010 (20 splits); S12 peak +0.163±0.036 at
  400×500, −0.0116±0.0109 at 1600×8000 (AF 0.863 there); s7 0.684/0.681/0.564/0.640/0.808
  (ctrl 0.504); s8 0.833/0.832 vs 0.845 (1dp 0.660/0.659 vs 0.622; PC 0.882); s1d
  0.976 vs 0.845 (+0.132); msprime 0.642/0.649; boundary table; s3b 3→0 probes.
- **Design integrity**: 3-way disjoint splits (sim exchangeable by construction; real
  panel shuffled within super-population before slicing); S10 orientation calibrated on
  the train half only (verified in code); CV meta-attacks out-of-fold; study-fitted
  attacks consume rel["loadings"] (the released object).
- **Honest scoping**: abstract and contribution 1 scope PC dominance to the simulated
  panels and state the real-data attenuation (0.652 < 0.688); Limitations repeats it.

## 1. Score: **7/10** for a top venue

## 2. Weaknesses (ranked) + 3. minimum fixes

**W1 (must-fix, text): the paper's real-data explanation is refuted by its own data.**
"Real allele-frequency spectra are less structured than the simulated panels" appears
twice (§7 real-data paragraph; §8 PC-channel paragraph). Measured directly from the
repo's panel: top-10 study-PC variance fraction **real 0.262±0.001 vs sim 0.144±0.003**
(matched 450×1363; 0.122 at the headline 900×3814 op point) and realized between-group
Fst-like statistic **real 0.135 vs sim 0.075**. Real is MORE structured on both metrics,
yet the PC channel attenuates there — the asserted explanation is backwards. Any
reviewer with data/external can check this in ten lines.
**Fix**: delete/replace both sentences with "the attenuation mechanism is not fully
explained" — or measure the plausible cause (real top-10 PCs dominated by super-pop
structure shared by members and non-members, so the in-sample-ness share of the
released spectrum is smaller).

**W2 (must-fix, text): truncated sentence.** Line 354 ends "…1\,dp/top-k---holds on"
and line 355 opens \paragraph{Real-data PC channel.}; the run-in heading makes the PDF
read "…—holds on Real-data PC channel. The released…" — an incomplete sentence missing
its object and period. **Fix**: complete it ("holds on real data as well.").

**W3 (one sentence + doc fix): a measured result is unreported, and misread in
CLAIMS_FROM_RESULTS.md.** s9a (0.7614±0.0116, s11_robust.json) appears nowhere in the
paper while the Discussion still hedges "a real batch-correlated QC table could
differ" — selective-omission exposure. Moreover the s9a design cannot measure what its
name claims: batch assignment is randomized independently of membership, so the
"which-SNPs-present" channel is structurally zero by design. Control run (this review):
uniform random 70% coverage gives AF AUROC 0.808/0.822; the 70/30 mixture gives
0.762/0.768 — i.e. 0.761 is the AF channel under heterogeneous (degraded) coverage,
not an added channel. **Fix**: one Discussion sentence reporting the measurement with
the correct interpretation; correct the claims-doc "new channel filled" wording.

**W4 (free strengthening): the real-data paired test costs nothing.** The claimed
"small consistent increment" (AF+PC vs AF on real chr22) is +0.0428, positive 20/20
splits, paired t p=1.3e-7, Wilcoxon/sign p=1.9e-6, bootstrap CI [+0.033, +0.053] —
computable from the existing s9 per-seed JSON. **Fix**: add the one sentence; closes
the claims-doc "real-data paired significance" missing-evidence item and strengthens
the real-data story (formally significant, if smaller, increment).

**W5 (minor text)**: (a) seed inventory (5/3/2) omits the new experiments (S9: 20
splits; S10: 25 seeds; S12: 12/cell; S13/S14: 20 splits); (b) "mid-marker corner" for
500 markers — the sparsest cell in the grid, not mid; (c) Fig 5(a) plots pc_naive under
both "PC naive" and "PC calibrated" labels (identical here only because orientation is
unanimous 20/20 — annotate or plot the train-calibrated value).

**Residual (disclosed, not blocking)**: one real dataset (chr22) only; PC channel not
cross-checked under msprime (S4 covers the AF tier only); haplotype channel future work.

## 4. READY? **Almost**

Yes after W1+W2 (must-fix text, ~15 min), W3 (one sentence + claims-doc fix), W4 (one
sentence), W5 (cosmetic). **No new experiments required — every fix is text-level.**
Do not submit with W1 in place: a reviewer with the artifacts can refute a factual
explanation in the paper, which would damage trust in the (otherwise verifying) numbers.
