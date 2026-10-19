## Round 1 (run_20261006_f68b0cd6) — Score: 7.5/10

### Raw Reviewer Response (verbatim)
Review complete. Full verdict written to `refine-logs/pi_review_verdict_round11_loop1.md`.

## What I verified (adversarially, artifact-by-artifact)
Read the paper, all six round-11 scripts, all six round-11 JSONs (plus shards/caches/maps), the four refine-logs docs, and the audit — then re-derived every round-11 number from its JSON key. The experimental substance holds: s39 enforcement fracs 3/3, governed 0.624/0.509; s40 bundle 0.962/0.947/+0.015 5-of-5, LOO 0.076, noise control 0.839; s41 planted 0.715/0.729/0.702, real 0.456/0.479/0.424 vs AF 0.866/0.864/0.838; s42 rho-sweep; s37 gap −0.03…−0.10; s38 A* 0.61, eps 23.9. Stale "remains argued, not built" removed; 15pp; ~20 round-8/9/10 regression spot-checks match.

## 1. Score: 7.5 / 10 (8.5 attainable after fixes)

## 2. Remaining weaknesses (ranked)
- W1 (Moderate): §governance "A* exceeds the binormal value 2AUC−1" is inverted; JSON gap_vs_2aucm1 negative on all five leaking channels, positive only on controls.
- W2 (Moderate): eps_lb≈24 clip-driven, not a sound finite-sample lower bound (b2 branch at TPR=1 with 1−TPR clipped at 1e-12); n_neg is 900, not ~300.
- W3 (Minor-mod): §replication says phased-IBD "remains future work" while §scale reports it run.
- W4 (Minor): Limitations points the biobank cohort at sec:governance; it lives in sec:scale.
- W5 (Minor): "base missingness shifts the absolute AUROC only weakly" — sweep is saturated, cannot characterize a weak shift.
- W6 (Minor): "full marker density … (50k markers strided uniformly)" self-contradictory.
- W7 (Minor, docs): ROUND11_DESIGN.md stale grid not refreshed; results-doc cohort CV 0.32 should be 1.32.
- W8 (Nit): s38 dead import anc_only_attack; governance.py docstring "remains argued".

## 3. Minimum fixes
W1 one-sentence swap; W2 state two-branch formula + n_neg 900 + ordering-only read; W3 delete/redirect future-work sentence; W4/W6/W8 one-line edits; W5 reword to "saturates regardless of base missingness"; W7 apply doc refresh + fix CV. No new experiments required.

## 4. READY? Almost — not as-is, yes after the fix list.

### Memory Update
- Verified: all six round-11 artifact groups match the paper; stale sole-egress sentence removed. Unresolved: paper-vs-JSON direction slip (A*), eps_lb estimator mismatch, §replication future-work contradiction.
- Pattern: round-11 additions to §scale inherited the plan's success criterion verbatim where it disagreed with the run's JSON; check direction of every "exceeds/positive" claim against the key.

---

## Round 2 (run_20261006_f68b0cd6) — Score: 8.0/10