# Round-5 fresh review — Round 1/4 (governance-layer round)

Reviewer: pi (senior privacy/ML reviewer persona). Method: direct artifact review —
main.tex read in full; governance.py, s21–s26, common.py read in full; every
checked number script-matched against `results/raw/*.json`; fresh compile
(12 pp, 0 errors, 0 undefined refs); Dirichlet-concentration behavior verified
empirically.

## What holds up (verified, not trusted)

- **All five previously-open limitations are addressed with real artifacts.**
  s21 Fst sweep (18 cells × 15 seeds, JSON verified), s22 admixture (2 arms),
  governance.py + s23 (real enforcement logic: undeclared/recipient/condition
  violations raise, over-precision coarsened, chain verifies; governed release
  0.622/0.510 = declared-safe tier), s24 (genuine 5-seed composition measurement,
  out-of-fold CV: cumulative 0.941 > max single 0.882, union > max on 5/5 seeds),
  s25 (pre-registered falsifiable prediction + null-family control:
  Δ_E=+0.2437 CI [0.2308, 0.2567] at p_san=1; null_family_paradox=False),
  s26 (planted-IBD positive control 0.866 vs unplanted 0.486; real chr22 0.463
  ≤ AF 0.653).
- **~40 paper numbers script-matched the raw JSONs this round** (s23 matrix and
  aggregates, s24 cumulative/max_single/flags, s25 paired CIs, s26 all cells,
  s21 grid and deltas, s22 deltas, abstract 0.88/0.84/0.94, per-seed combined>AF
  on 25/25 seeds with combined min 0.9217 > AF max 0.8639).
- Regression check: popets2026 → Berrang et al.(2026) ✓; fedlmm2025 authors ✓;
  "canonical" scoping ✓; msprime 0.642/0.649 honest reporting ✓; s19/s20
  untouched (mtimes) ✓; no `??`; no TODOs.
- Scope compliance: prototype honestly disclaimed (plaintext, release-object
  layer; sole-egress/SPU binding argued, not built) in Limitations and in code
  docstrings. No network/remote-exec code in the new scripts.

## 1. Score: 8/10

## 2. Critical weaknesses (ranked) + 3. minimum fixes

**W1 (must fix — interpretation inverted).** The s22 "discrete structure"
conclusion is backwards relative to its own experiment. Paper (abstract,
contribution 5, §7, Limitations): "under continuous admixture (Dirichlet
weights, α=0.2) the PC channel falls below AF (−0.115), recovering dominance
only as the structure becomes discrete (α=10: +0.096)"; "the effect requires
discrete structure". For a symmetric Dirichlet, α=0.2 gives *sparse, near-one-hot*
weights (verified: mean modal fraction 0.81) — i.e., the *most discrete-like*
arm — and α=10 gives *near-uniform* weights (modal fraction 0.43) — the *most
admixed* arm. So the data say: PC dominance is recovered in the MOST ADMIXED
setting and lost in the most discrete-like one — the opposite of the claim as
written. The result actually aligns with the s21 structure-strength story (gap
grows as true structure weakens), not with a discrete-vs-continuous story; and
the design cannot isolate "discrete vs continuous" anyway (both arms are
continuous-weight mixtures differing only in concentration).
*Minimum fix:* re-interpret s22 in text — "as admixture weights become more
uniform (structure weakens), the PC−AF gap grows, consistent with the Fst
sweep; this experiment varies admixture concentration and does not isolate
discrete-vs-continuous structure" — and strike "needs/requires discrete
structure" from abstract, contribution 5, §7, and Limitations. (If the
discrete-structure claim is wanted, it needs a new experiment: α→0 one-hot arm
vs a continuous arm at matched structure strength.)

**W2 (must fix — nominal/measured conflation).** §7: the Fst sweep "crosses
zero near F_ST≈0.10 … the structure level of the real 1000G panels" (and Fig 7a
shades 0.09–0.15 on the nominal axis). The x-axis is the *nominal* simulation
parameter; the real panels' ≈0.10 is the *measured* Fst-like statistic. By the
paper's own numbers (and s21's fst_like: measured ≈0.052 at nominal 0.08), the
measured-structure crossover is ≈0.04–0.06 (n1600/n800), i.e. BELOW the real
0.10–0.14 — the real panels sit *beyond* the crossover. This actually
strengthens the attenuation story; the claimed coincidence is an
apples-to-oranges alignment. Also "decreases monotonically at every site size"
is only strictly true at n800/n1600; the n300 row wiggles (+0.0061→+0.0074,
within the code's 0.01 tolerance) and is at ceiling (both channels ≈0.99).
*Minimum fix:* report the crossover in measured fst_like terms (≈0.04–0.06) and
state the real level lies beyond it; qualify n300 as ceiling-bound; move/label
the Fig 7a shade accordingly.

**W3 (must fix — claim vs code).** The "phased-haplotype/IBD attack" is not
phased. real_phased.py keeps H0/H1 on disk, but s26's real arm sums them back
to dosage (`H0[si]+H1[si]`) — the code comment "and to phased for the real one"
is unfulfilled; neither arm uses phase. Also the mechanism is a count of
exactly-matching 8-marker windows, not "long shared runs of identity-by-state".
The honest conclusion (0.99 = fitted-table artifact; real sharing among
unrelated 1000G samples leaks little) survives — the planted control and
real-data chance result are real — but the label overstates.
*Minimum fix:* reword to "identity-by-state window-sharing attack on unphased
dosages (the phased VCFs are loaded but phase is not yet exploited)", or
actually run the attack on H0/H1 haplotypes. (s26's docstring also promises a
15kb-thinning sensitivity control that is not implemented — fix the docstring.)

**W4 (should fix — self-fulfilling flagging).** The composition audit's
"anti-composition clause flags the refinement" is tautological as implemented:
s24 writes `meta={"refines": "af"}` into the ledger records itself, then
`Ledger.relations()` reads that metadata back and "flags" it. The monitor does
not enforce anti-composition; the flags are producer-asserted, not inferred
from released content. *Minimum fix:* one sentence in §7 ("the relation flags
are producer-asserted metadata, not inferred from content; content-based
relation detection is future work") — or make the check content-based.

**W5 (should fix — novelty honesty + missing control side).** The s25 paradox
result is a corollary of the already-measured exposure curve: Δ_E = e(p_stop) −
e(p_san); the seven analyst families only determine the stopping rung, and
quality parity is trivial (the quality metric barely varies with k, so
"equal GWAS quality" holds by construction). Unreported: for the null family at
p_san≥2, Δ_E = −0.24/−0.25 — a sanitized block *finer* than the analyst's need
*raises* exposure (the sanitized release is itself the leak). *Minimum fix:*
one sentence acknowledging the result restates the precision frontier through
the analyst-decision wrapper, plus one sentence reporting the negative-side
control.

**W6 (minor).** Seed inventory: s21 (15 seeds, given inline), s22 (10 seeds
parametric, 3 msprime), s26 (5 seeds sim, 10 real splits) are not in the
§5 operating-point seed list. "Loadings≫eigenvalues holds at every point" of
the s21 grid — the sweep has no eigenvalue-only arm (it exists at the s22 and
standard operating points). Tracker says "10 pages"; fresh compile is 12.
*Minimum fix:* add the seed counts; drop "at every point" or scope it; sync the
tracker page count.

## 4. READY? Almost

All headline and governance numbers verify against artifacts; all five prior
limitations are closed with real experiments. What blocks "ready" is textual
claim accuracy, not data: W1 (inverted interpretation, in the abstract), W2
(nominal/measured conflation), W3 ("phased" label). All three are paragraph-
level rewrites; no re-runs needed except optionally W1's follow-up experiment.
After W1–W3 (and ideally W4–W6), this is submission-ready.
