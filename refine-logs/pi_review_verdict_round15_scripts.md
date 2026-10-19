# Round-15 limitation-closure scripts — code review verdict

**Scope**: s51_real_lowstruct.py, s52_haplo_second.py, s53_eigbasis.py, s38_privacy_metric.py (rule-of-three edit), s37_fitsens.py, against EXPERIMENT_PLAN_LIMITS.md; s48 helpers reused by them. All claims below verified by reading the code AND by empirical probes on the cached data (planted-control power curve, label alignment, LSR saturation, within-pool Fst).

**Verdict**: s51 correct as-is (two wording-level minors). s38's eps machinery and floor pairing are correct (two minors). s37 correct (one minor). s53's bases are fixed and the orientation is held-out, but it does not implement what B4 planned (MAJOR, claim-level). s52's real arm is correct, but its **planted positive control does not recover — CRITICAL** — and the plan's B3 premise ("planted-positive control still recovering", cited from s36) is stale: the current s36/s32 LSR implementation union-saturates and cannot reproduce the cited numbers (MAJOR, adjacent).

---

## CRITICAL

### C1. s52 planted positive control fails (planted AUC 0.525 chr22 / 0.536 chr1 ≈ chance)
The script's own success criterion ("at least one detects the planted segment") and B3's ("the planted-positive control still recovering") are **not met by the recorded results**. Verified root cause: **statistical power, not a coding slip in the self-exclusion**.

- Same-ancestry pairwise IBS share at this cache (m=4863, MAF≥0.05, 1.5 kb thinning) is ≈0.70–0.77; the max over ~478 released haps sits at ≈0.789.
- A 300-marker exact segment lifts a carrier's share by only `seg·(1−p_match)/m ≈ 300·0.23/4863 ≈ +0.014` — **inside the max-noise floor**. Measured: planted-member mean 0.7915 vs unplanted 0.7874 (carrier-vs-noncarrier AUC within members: 0.549).
- Design wrinkle that compounds it: `plant_phased_ibd(members_frac=0.6)` modifies **60% of the whole panel** (both future members and non-members, contiguous in panel order), but non-members are then scored on the **original H0** while members are scored on P0. This accidentally produces the intended class asymmetry (members 60% carriers, non-members 0%), but under the plant's own ground truth the non-member carriers' true haplotypes are P0[ni]; scoring them on H0 is inconsistent.
- **Fix is not sufficient by si-only planting alone** — measured: si-only planting at seg=300 still gives AUC 0.521. Power curve (si-only planting, members only carry, same statistic/split): seg=300 → 0.52, seg=1000 → 0.70, seg=1500 → 0.80, seg=2000 → 0.80.

**Exact fix** (both parts required):
1. Plant only into the split's members, after the split: `Q0 = H0.copy()`; per super-pop `mem = np.intersect1d(si, np.where(sup==p)[0])`, donor = `mem[0]`, recipients = random 60% of `mem[1:]`; copy `H0[donor, st:st+seg]`. Then `P0==H0` on all non-member rows, the ground-truth inconsistency disappears, and `snP = score_targets(ni, H0, H1, relQ, self_non)` needs no change.
2. Scale the segment to the panel: `seg = max(300, H0.shape[1] // 5)` (≈970 for m=4863; 1500 gives AUC 0.80). Report the power bound honestly: the whole-hap IBS-share statistic cannot certify detection of segments below ~20% of the panel; short-segment claims must stay scoped to the run-based estimator.
3. With seg≈1500–2000 the `n_close` count also fires on carriers (share 0.85 > 0.78), so the fix revives both statistics; rerun both chromosomes.

## MAJOR

### M1. s53 does not implement B4's "eigenvalue-informed" basis (claim-level)
The identity/standard-axis and DCT-II bases are genuinely **fixed and data-independent**, and both land at chance (0.49–0.52 vs full 0.81–0.88) — but they carry **no eigenvalue information** (the script's own docstring concedes "not tied to the genotype covariance"). What is closed: "the negative is not an artifact of basis randomness." What is NOT closed: L4's planned probe of an *eigenvalue-dependent* attack. Note the plan's specific construction is ill-posed (eigenvectors are unidentifiable from eigenvalues + isotropy — any orthogonal basis fits), which is presumably why the substitution was made.
**Exact fix** (choose one, and make the paper text match):
- (a) Word the paper as "fixed, data-independent orthonormal bases (standard-axis, DCT) also fail — the failure does not depend on the basis being random", and drop the words "eigenvalue-informed"; or
- (b) add the actual eigenvalue-dependent statistic, one line over the existing machinery: score `d(x) = Σ_j ((x−μ)·b_j)² / λ_j` with **released** top-k eigenvalues and the fixed basis — an attacker who knows only the spectrum and uses fixed directions. If that is also at chance, L4 closes as planned.

### M2. (adjacent) The plan's B3 baseline numbers are stale: current `lsr_score` union-saturates
The cited s36 results (real-arm LSR ≈343, planted AUC 0.671) **cannot be produced by the current code**: the s36 JSON's schema (`agg.real/planted.{auc,mem,non}`) does not match the current `s36_phased_selfexcl.py` output keys (`lsr_auc/af/...`, `planted/unplanted`), and the current `s32_phased_ibd.lsr_score` computes the per-marker **any-hap union** before run-length: measured on the s36 cache, an *excluded unrelated* target vs only 120 released haps scores **4088/4431** (saturates), while the true per-single-hap longest shared run is 389. A re-run as-is would report ≈m for everyone and a planted control ≈0.5 for the wrong reason.
**Exact fix**: in `lsr_score`, take runs **per released haplotype** then max over haps — replace `_runs_maxlen(eq.any(2))` with a per-column `max(_runs_maxlen(eq[:, :, j]))` over the block axis (keep the 50-hap chunking) — then re-run s36's positive control + real arm before the paper cites "the existing LSR estimator" anywhere, and only then pair it with s52. (If the OR-union was intentional as a set-cover statistic, it must be renamed and its saturation documented — but it cannot serve as the "longest-shared-run" estimator the paper describes.)

## MINOR

### m1. s51 pooled arm: reported `fst_like` = 0 by construction, not a measurement
Single-group label ⇒ `fst_like` returns exactly 0. The hidden within-pool structure is real but small (measured pop-level Fst on the s51 split: EUR 0.016, EAS 0.009, AFR ≈0/−0.019 noise) — genuinely low-structure vs the ladder's 0.088, but a table showing "0.0" invites the "you labeled the structure away" objection. **Fix**: `_sup_label` already loads `lab` (pops) and discards it; report `s48.fst_like(site, lab[si])` alongside for the pooled arm (one line in the Arm A loop).

### m2. s51 `POOL_PER=160` is below the per-group cap the docstring claims
Caps are 167 (EUR), 168 (EAS), 220 (AFR). Docstring says "at the largest per-group n the pool allows". **Fix**: either `per = G5[m].shape[0] // 3` per super-population, or reword to "a uniform per-group cap of 160 (panel caps 167/168/220)". The panel bound itself is honestly reported via `n_pool`.

### m3. s52 `n_close` at THR=0.78 is vacuous at this panel length
share > 0.78 requires ≈3790/4863 matching markers — unreachable except cryptic relatives; the statistic is ≡0 almost everywhere, so its "near chance" AUROC is tautological, and it misses the plant too (planted n_close 0.526). **Fix**: drop the statistic, or scale THR with the plant/panel (with the C1 fix's seg≈1500 it becomes non-vacuous as-is); do not quote its real-arm 0.496 as evidence of anything.

### m4. s52 docstring promises a statistic the code doesn't compute
"A max-over-window statistic is also reported but is saturated ... documented as such" — no such statistic exists in `run()` or the JSON. **Fix**: delete the sentence (the saturation point is already made by M2's evidence) or actually compute it.

### m5. s52 vs plan: "both on the same members/non-members" is not implemented
The new estimator uses `per=80`, `ld_window=1500, max_snps=30000`; the LSR arm used `per=150`, `ld_window=1000, max_snps=20000` — different splits and caches, so the appendix "per-estimator" comparison is unpaired. **Fix**: reuse `s32_phased_ibd.split_idx` and the same phased cache in s52 (after M2's fix), or present the rows as separate experiments, not a paired comparison.

### m6. s38 `eps_order_inversion` note is hardcoded prose
It matches the current numbers (af_full 0.8446 / 4.12 vs af_2dp 0.8396 / 4.19 at δ=0.01; 3.41 vs 3.55 at δ=0.05 — verified), but a rerun that flips the ordering would silently leave a false claim in the results JSON. Also ΔAUC = 0.005 is within seed noise at n=900/900 — the example is a metric-property illustration, not a measured effect. **Fix**: compute the direction (`a["auroc"] > b["auroc"] and a["eps_lb_0.01"] < b["eps_lb_0.01"]`) before writing the note, and store per-seed stds for the pair.

### m7. s38 `valid` is now consistency-by-construction (expected, but don't oversell)
With matched floors the flag can only fail numerically — it certifies that the reported ε satisfies the inequalities on the same grid, not an independent check. The docstring's framing is already correct; keep the paper's wording descriptive ("internally consistent at the stated resolution floors"), never "verified DP".

### m8. s37 `full` rung omitted at (16000, 3000) contrary to its own feasibility rule
The docstring's rule (n_site²·4B ≤ ~20 GB) allows n_site=48000 (9.2 GB float32 Gram), and the plan's sweep is {5000, 10000, 20000, full}. First-cell stability is already good (gap −0.0300 → −0.0296 across 21%→100% fit fraction, sd 0.0078); the (16000,3000) sweep was still mid-run at review time. **Fix**: run the nfit48000 rung if the GPU allows, else state in the docstring/plan that `full` was not run at 48000.

---

## Verified correct (no action)

- **s51 Arm B ladder (Q1)**: composition genuinely held fixed — every rung is exactly 1:1:1 AFR:EAS:EUR (panel thirds 220/167/168 ≥ 160, so `q=min(per, size//3)=per` at every rung); measured `fst_like` stable across rungs (0.091→0.087); reproduces the plan's probe numbers exactly (+0.120/+0.042/+0.007/−0.026); chr1 also monotone-crossing (+0.058/+0.002/−0.010/−0.018). Disjoint 3-way split per group; controls at chance in every rung (refb 0.47–0.51, eigo 0.50–0.53).
- **s51 Arm A (Q1)**: a genuinely low-structure real cell (hidden pop-level Fst 0.009–0.016, the lowest 1000G supports), correct disjoint split, controls at chance. PC>AF persists at n=160 as the plan predicts.
- **Label alignment (s51/s53/s48)**: `real_multi._parse_vcf` emits G rows in panel-file filtered order and nothing downstream reorders; `load_for`/`_sup_label`/`crit` all derive labels from that same order. Verified in the loader source; no misalignment.
- **s52 real arm + leave-self-out (Q2)**: self-rows `[t, n_st+t]` match the `vstack([H0[si], H1[si]])` layout exactly; `setdiff1d` removal is correct; non-members correctly have no self-rows (their haplotypes are not in the release); the middle third is unused, so members/non-members are disjoint. Real-arm negative valid: 0.488±0.011 (chr22), 0.497 (chr1). This part of B3 stands.
- **s53 fixedness + orientation (Q3)**: both bases are deterministic and data-independent (`_identity_basis`, `_dct_basis` — no input dependence); orientation is fixed on the held-out calibration half only via the shared `cal` mask (seed+100), the same mask for all four bases, so the comparison is fair and no scored row chose its own sign. Results: identity 0.49–0.52, DCT 0.50–0.51 vs full 0.81–0.88 — at chance, matching the random-basis ablation.
- **s38 eps machinery (Q4)**: `eps_lower_bound`'s two branches are exactly the (ε,δ)-DP ROC inequalities (TPR ≤ e^ε·FPR+δ and its complementary-direction twin 1−FPR ≤ e^ε(1−TPR)+δ; the remaining two are covered because the evaluated grid is closed under threshold complement); linear-fractional monotonicity means grid endpoints suffice between points; `validate_tradeoff`'s `fnr ≥ g` is algebraically the same two constraints. **The floors match at every call site**: `eps_lb_*_valid` uses 1/n_pos in both the bound and the check; `eps_lb_rule3_*_valid` uses 3/n_pos in both; same `alpha_floor=1/n_neg` and same `sel` grid. The rule-of-three floor (3/n_pos, zero observed FNs in n_pos trials) is the correct confidence bound and floors in the conservative direction (shrinks b2, lowers the claimed ε). Results confirm: all channels valid=True including rule3, top channel af_sitewise ε≈5.65 at δ=0.01 — B2's success criterion met.
- **s37 (Q5)**: `n_fit` semantics are exactly "PCA fit on a random subset of the study rows" (`fit_idx` drawn from site only; AF channel unaffected); SWEEP keys (8000,8000)/(16000,3000) are the two largest delivered cells in `s37_realsynth_scale.json`; the delivered `n_fit`/`fit_frac` bookkeeping is consistent.
- **s48 helpers as reused**: `_orient_auroc` (calibration-half sign, held-out evaluation, shared mask across channels) and `cohort_cell`'s 3-way disjoint split behave as the established discipline requires; no self-inclusion is introduced beyond the intended released-object semantics (study-fitted loadings/AF contain members by construction — the paper's documented channel).

## Bottom line
- Merge-ready as measurement code: s51, s38, s37.
- s53: merge after the paper text is scoped to "fixed data-independent bases" (or add the eigenvalue-weighted variant, M1b).
- s52: **do not write the B3 appendix row until C1 is fixed and rerun** — the real-arm negative stands, but "the planted-positive control still recovering" is currently false in the recorded output, and the LSR baseline it would be compared against (M2) is itself unreproducible with current code.
