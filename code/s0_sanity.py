"""S0 sanity: verify the pipeline + attack meter before trusting any leakage
number. Paper-validation mode (CPU only).

Checks:
  1. population stratification inflates uncorrected lambda_GC; PCA correction tames it.
  2. a planted signal is detectable.
  3. the AF-LRT attack beats chance on the F5 endpoint, with the site vs a
     population-matched, disjoint reference set.
  4. the ancestry-only control attack is weaker than the real attack.
  5. the attack is ~chance on random-noise input (meter sanity).
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    simulate_population, split_population3, release_objects,
    marginal_association, genomic_control_lambda, allele_freq,
    af_lrt_attack, anc_only_attack, auroc, tpr_at_fpr, pca,
)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")
os.makedirs(OUT, exist_ok=True)


def run(n_site_per_group=300, n_ref_per_group=300, n_groups=3, n_snps=4000, seed=0):
    pop = simulate_population(n_groups=n_groups, n_per_group=1000, n_snps=n_snps, seed=seed)
    # three-way split: study members, attacker reference, independent non-members
    site, ref, nonmem = split_population3(pop, n_site_per_group, n_ref_per_group,
                                          n_ref_per_group)
    rel = release_objects(site)

    # --- 1/2 utility sanity ---
    mask = rel["mask"]
    _, p_uncorr = marginal_association(site["G"], site["case"], n_pcs=0, correct_pcs=False)
    _, p_corr = marginal_association(site["G"], site["case"], n_pcs=10, correct_pcs=True)
    lam_uncorrected = genomic_control_lambda(p_uncorr)
    lam_corrected = genomic_control_lambda(p_corr)
    causal = site["causal"]
    mask = rel["mask"]
    causal_in_mask = [c for c in causal if mask[c]]
    min_p_causal = (float(np.nanmin(rel["pvals"][causal_in_mask]))
                    if causal_in_mask else float("nan"))

    # --- 3/4 attack sanity (F5 endpoint = allele frequencies) ---
    members = site["G"]
    nomembers = nonmem["G"]
    study_af = allele_freq(members)
    ref_af = allele_freq(ref["G"])
    scores_m = af_lrt_attack(members, study_af, ref_af, mask=mask)
    scores_nm = af_lrt_attack(nomembers, study_af, ref_af, mask=mask)
    scores = np.concatenate([scores_m, scores_nm])
    labels = np.concatenate([np.ones(len(scores_m)), np.zeros(len(scores_nm))])
    a = auroc(scores, labels)

    # structure-only control: project BOTH target panels onto the study-fitted
    # PC basis, distance to the study centroid. Members (in-sample w.r.t. the
    # study) sit farther; this is the pure structure/ancestry channel.
    _, _, load_s = pca(site["G"], k=10)
    mean_s = site["G"][:, mask].mean(0)
    anc_scores = np.concatenate([
        anc_only_attack(members, load_s, mean_s, mask),
        anc_only_attack(nomembers, load_s, mean_s, mask),
    ])
    a_anc = auroc(anc_scores, labels)

    # noise meter
    rng = np.random.default_rng(seed + 7)
    Gm = rng.integers(0, 3, size=members.shape).astype(np.int8)
    Gn = rng.integers(0, 3, size=nomembers.shape).astype(np.int8)
    s_m = af_lrt_attack(Gm, study_af, ref_af, mask=mask)
    s_n = af_lrt_attack(Gn, study_af, ref_af, mask=mask)
    a_noise = auroc(np.concatenate([s_m, s_n]),
                    np.concatenate([np.ones(len(s_m)), np.zeros(len(s_n))]))

    res = dict(
        config=dict(n_site=members.shape[0], n_ref=nomembers.shape[0],
                    n_snps_tested=int(mask.sum()), seed=seed),
        lambda_gc_uncorrected=float(lam_uncorrected),
        lambda_gc_corrected=float(lam_corrected),
        n_causal=int(len(causal)), n_causal_tested=int(len(causal_in_mask)),
        min_p_causal=min_p_causal,
        attack_auroc=float(a),
        ancestry_only_auroc=float(a_anc),
        noise_auroc=float(a_noise),
    )
    print(json.dumps(res, indent=2))
    with open(os.path.join(OUT, "s0_sanity.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
