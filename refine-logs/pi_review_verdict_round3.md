# Cross-model review verdict — round 3, paper-level (2026-10-02)

Reviewer scope: PoPETs/USENIX. I read `paper/main.tex`, verified every number
in the paper against `results/raw/*.json`, and verified in code+data that all
five round-2 must-dos were actually implemented this time (S3b, S4 re-run,
S1c, S6, S2c). All Table 1 numbers match s1c exactly (including the
per-ancestor caption 0.743/0.746/0.709); the mechanism claims match s2c; the
boundary table matches s5; the dimension-sweep sentence matches s6.

## Round-2 must-dos: verified ledger

| must-do | verified? | evidence |
|---|---|---|
| S3 implemented for real | **YES** | `s3b_paradox.py`: stopping rules computed from the revealed spectrum (scree threshold / largest gap), quality measured by re-running the association pipeline (calibration + recovery), probes emergent (3→0 both families, 3 seeds), disjoint negatives. Residual assumption, admitted in the docstring: the AF-precision→resolvable-PCs mapping is a modeling choice; quality parity is a consequence of the stability-stopping design (stability ⇒ same k ⇒ same measured quality), and the measured exposure reduction is only ~0.01 AUROC (0.857→0.845) because 2dp-vs-3dp AF barely differ. The paper's "conjecture, not a result" framing is honest and accurate. |
| S4 re-run with disjoint split | **YES** | `s4_msprime.py` now scores a disjoint non-member panel; 0.642 full / 0.649 1dp (was 0.716/0.627 with the biased set) |
| S1c 5 seeds + per-ancestor | **YES** | `s1c_frontier.json`: 5 seeds, per-ancestor ≈uniform (0.709–0.746 at F5) |
| S6 dimension sweep | **YES** | `s6_dim.json`: 9 operating points, Δ always positive (0.031–0.124); **pcdist-only beats the AF base at all 9** — this is the strongest support for the headline |
| S2c mechanism formalised | **YES** | matches my round-2 diagnostic: pcdist-only 0.791, pcalign 0.504 (chance), T2 Δ +0.003±0.006, T3 Δ +0.103±0.010, eig ratio 0.9998 |

---

## 1. Is the headline claim ("diagnostics release dominates") supported?

**SUPPORTED**, with one required rewording.
The dominance direction is real and now robustly evidenced: the PC-projection
attack alone (0.79) beats both AF summary releases (0.75 at 2 dp, 0.72 at full
precision) and does so at **all nine** dimension-sweep operating points
(pcdist-only > base in every row of s6); combined they reach 0.86; per-site AF
(0.858) is a still-larger single channel; the mechanism (in-sample-ness of the
study-fitted spectrum) is identified with the ref-fit control collapsing the
increment to ≈0 (+0.003±0.006 over 5 seeds) while attacker-realistic centering
leaves it intact (+0.103±0.010).

**But one abstract clause is falsified by the paper's own numbers:** (ii) says
the PC diagnostic is "the single largest leakage channel" while (iv) correctly
reports per-site AF "worse still (0.86)" — and 0.858 > 0.791. The PC spectrum
is the largest *intermediate diagnostic*; per-site AF is the largest single
channel. Reword (ii) accordingly; as written, a careful reviewer flags the
self-contradiction in the abstract.

Two smaller precision fixes: (a) "full allele-frequency summary release (0.75)"
mislabeled — 0.753 is the **2 dp** release; full precision is 0.720 (fixing
this makes the PC claim stronger: 0.79 > 0.75 > 0.72); (b) "magnitude scales
with the SNP-to-member ratio" — the sweep shows Δ depends on site size too
(ratio 30 → +0.072 vs ratio 20 → +0.120); say "ranges from +0.03 to +0.12
across operating points, largest at small sites, always positive."

---

## 2. Self-reference-bias correction: real contribution or self-inflicted-repair?

**Both — real but minor, and self-inflicted in origin.** The quantified
demonstration (≈0.2 AUROC inflation; 0.92→0.72) of a specific, easily-committed
control-design error — scoring the attacker's own reference individuals as
negatives, which turns the one-sided membership test into an easier two-panel
distinguishability test — is a genuine methodological caution. The MI
literature warns about ancestry confounds but rarely quantifies this particular
bias. It strengthens the paper's credibility as honest accounting. It should
NOT be contribution #1 (currently is): a PoPETs reviewer will read it as due
diligence, and it originated as this project's own round-1 flaw. Demote it to a
methodology subsection; it cannot carry a paper, but documented alongside the
corrected measurements it adds real value.

---

## 3. Single biggest remaining weakness

**All-simulated: no real cohort and no real published summary-stat file.** The
internal machinery is now sound, honest, and internally consistent — disjoint
controls, mechanism identified, boundary curve, dimension sweep, emergent-analyst
conjecture properly framed. What remains is external validity: every magnitude
is simulator-specific (the paper says so plainly), and the governance
implication ("treat the PC spectrum as a governed object") rests entirely on
simulation. This is not fixable within the current evidence and is the gap a
PoPETs reviewer will weigh heaviest.

Most dangerous *fixable* flaw, named separately because it is a trust issue:
**the msprime cross-check is cited in the Data paragraph but its numbers are
nowhere in the paper — and they now complicate the story.** The re-run shows
full 0.642 vs 1dp 0.649: the 1dp rounding benefit **vanishes** at the msprime
operating point (within noise), while the boundary section claims "leakage
collapses with ... coarse (1 dp) rounding" from the BN simulator alone.
Omitting a cited cross-check whose numbers would qualify the narrative is
exactly the pattern reviewers punish hardest. Report it honestly — framed as
"at weakly-structured demographies the quantization effect is within noise,"
it actually *reinforces* the saturation/coarse-cutoff story.

---

## 4. Most likely rejection reason; fatal or fixable?

**Most likely rejection:** "All simulation; the qualitative phenomenon (AF
release → membership leakage) is 15+ years old (Homer, Sankararaman); the
novel quantitative claim's mechanism is in-sample-ness of a fitted PCA, whose
magnitude is an artifact of the simulator's dimensionality; the second attack
family is non-functional by construction; no real cohort or real summary-stat
file."

**Fixable within current evidence: partially.** The presentation-level clauses
are all fixable now (the abstract self-contradiction, the unreported msprime
cross-check, Table 2's caption mismatch, missing citations, thin threat model
— see the checklist below). The all-simulated clause is not fixable without
new work; the paper's defense is the one it already takes — the methodology
(the controlled protocol + mechanism identification) is the contribution, and
the limitations are stated plainly. That scoping is what makes the submission
defensible rather than the numbers themselves.

**Pre-submission checklist (all text edits + one tiny code addition; no new
modules):**
1. Reword abstract (ii): PC diagnostic = largest *intermediate diagnostic*
   channel; per-site AF is the overall ceiling. Fixes the (ii)/(iv)
   contradiction.
2. Report the msprime cross-check (0.64/0.65) in the results or boundary
   section, including the null rounding effect at that operating point.
3. Fix Table 2: caption says 5 seeds, numbers are S2b's 3-seed values
   (0.753/0.861/+0.108); S2c's 5-seed values are 0.758/0.860/+0.102. Use the
   5-seed numbers.
4. Fix the T3 range: paper says "+0.10 to +0.11"; per-seed is 0.083–0.112
   (say "+0.08 to +0.11, mean +0.10").
5. Add the label-permutation control to s2c (the paper claims "returns chance
   (≈0.50)"; it currently lives only in `results/raw/_s2b_diagnostic.py`,
   not the official pipeline), or cite that script.
6. Relabel "full allele-frequency summary release (0.75)" → "the AF summary
   release (0.75 at 2 dp; 0.72 at full precision)".
7. Reword the dimension-sweep sentence per §1(b).
8. Complete the related-work citations (the paper currently has none —
   `\todo{full citations}`), and tighten the threat model (A2 "constrained
   workflow author" and A3 "adaptive analyst" need one concrete sentence each).
9. Demote the self-reference correction from contribution #1 to a methodology
   subsection; contribution order should lead with the diagnostics finding.

---

## 5. Verdict: submission-plausible as a measurement paper?

**Yes — conditionally.** All five round-2 must-dos were genuinely implemented
(this time verified in code and data, not prose); the evidence base is sound,
honest, and internally consistent; the mechanism is identified with proper
controls; the limitations are stated without spin. With the nine-item
checklist above (a day of text edits and a ten-line code addition, no new
experiments), the paper is a legitimate measurement-methodology submission:
novel question (intermediate diagnostics' marginal leakage, with the PC
spectrum identified as the dominant intermediate channel via an
in-sample-ness mechanism), clean controls, honest negative results (λ_GC and
condition bit at exactly zero), and an actionable boundary curve.

For a PoPETs/USENIX main track it remains borderline on the all-simulated
clause — that is the one thing the current evidence cannot cure, and it should
be the next round of work (a real released AF vector, e.g. 1000 Genomes
populations as pseudo-sites) rather than more simulation. As a measurement
paper with the scoping it already has, it is submission-plausible now.
