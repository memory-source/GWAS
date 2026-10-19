"""S20 - two corrections/analyses flagged in review.

(1) M3 missingness as a genuine channel. The earlier batch-missingness variant
    (s11/s15) randomized batch assignment independent of membership, so its
    AUROC was coverage-degraded AF-LRT, not a missingness channel. Here M3 is
    implemented as the standard genomic-missingness channel: each individual has
    a called/missing pattern across markers; the study releases a per-SNP call
    rate table; a target's own called-set likelihood is scored under the study's
    released call-rate table vs the reference's. This measures leakage of the
    missingness table *as a released object*.

(2) Structure comparison and the Fst expectation. Under Balding-Nichols the
    between-group variance ratio Var_k(p_k)/(p_bar(1-p_bar)) equals the Fst
    parameter, so the simulated panel's nominal Fst is a direct target; the
    real 1000G super-population value is emergent. We report both so the
    "real panel is less/more structured than simulation" comparison is anchored
    to a theory value rather than a single statistic.

Output: results/raw/s20_extra.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, auroc)
from real_multi import get_real_panel_chr, three_way_split

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def call_pattern_channel(seed, missing_member=0.05, missing_non=0.03, holdout=True):
    """M3 channel, with the self-reference fixed. The released per-SNP call-rate
    table is estimated from a TRAIN half of the study; targets (the other half
    of members + the non-members) are scored out-of-sample. Members are genotyped
    in a study batch with a slightly higher per-marker call-failure rate
    (``missing_member`` vs ``missing_non``); the released call-rate table then
    differs between study and reference, and a target's own called-set pattern
    leaks membership. With missing_member == missing_non the table carries no
    membership signal (the null).
    """
    rng = np.random.default_rng(seed)
    pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
    site, ref, non = split_population3(pop, 300, 300, 300)
    m = site["G"].shape[1]

    def patterns(n_ind, rate):
        return rng.random((n_ind, m)) > rate

    cm = patterns(len(site["G"]), missing_member)
    cn = patterns(len(non["G"]), missing_non)   # non-members: reference-batch rate
    cr = patterns(len(ref["G"]), missing_non)
    if holdout:
        tr = rng.random(len(site["G"])) < 0.5
        study_call = np.concatenate([cm[tr], cr]).mean(0)   # table excludes test members
        members = cm[~tr]
    else:
        study_call = np.concatenate([cm, cr]).mean(0)
        members = cm
    ref_call = cr.mean(0)
    eps = 1e-3

    def ll(called, table):
        return (called * np.log(table + eps) + (~called) * np.log(1 - table + eps)).sum(1)

    d_m = ll(members, study_call) - ll(members, ref_call)
    d_n = ll(cn, study_call) - ll(cn, ref_call)
    y = np.concatenate([np.ones(len(members)), np.zeros(len(cn))])
    return float(auroc(np.concatenate([d_m, d_n]), y))


def fst_like(G, groups):
    af = G.mean(0) / 2.0
    hs = np.mean([2 * (G[groups == g].mean(0) / 2.0) * (1 - G[groups == g].mean(0) / 2.0)
                  for g in np.unique(groups)], axis=0)
    ht = 2 * af * (1 - af)
    return float(np.nanmean((ht - hs) / np.where(ht > 0, ht, np.nan)))


def m3_sweep(seeds=(0, 1, 2), base_rate=0.03, deltas=(0.0, 0.005, 0.01, 0.02, 0.05)):
    """Sweep the study-batch call-rate excess so the missingness channel's
    magnitude is interpretable rather than a single ~1.0 point. Reuses one
    simulated panel per seed (the channel is genotype-independent)."""
    out = {}
    for d in deltas:
        vals = []
        for s in seeds:
            rng = np.random.default_rng(s)
            pop = simulate_population(n_groups=3, n_per_group=600, n_snps=4000, seed=s)
            site, ref, non = split_population3(pop, 200, 200, 200)
            m = site["G"].shape[1]
            cm = rng.random((len(site["G"]), m)) > (base_rate + d)
            cn = rng.random((len(non["G"]), m)) > base_rate
            cr = rng.random((len(ref["G"]), m)) > base_rate
            tr = rng.random(len(site["G"])) < 0.5
            study_call = np.concatenate([cm[tr], cr]).mean(0)
            ref_call = cr.mean(0)
            eps = 1e-3
            ll = lambda c, t: (c * np.log(t + eps) + (~c) * np.log(1 - t + eps)).sum(1)
            members = cm[~tr]
            dm = ll(members, study_call) - ll(members, ref_call)
            dn = ll(cn, study_call) - ll(cn, ref_call)
            y = np.concatenate([np.ones(len(members)), np.zeros(len(cn))])
            vals.append(float(auroc(np.concatenate([dm, dn]), y)))
        out[f"delta_{d}"] = [float(np.mean(vals)), float(np.std(vals))]
    return out


def run():
    res = {}
    # (1) M3 channel
    print("=== M3 missingness channel ===", flush=True)
    m3 = [call_pattern_channel(s) for s in range(10)]
    res["m3_channel"] = dict(per_seed=m3, mean=float(np.mean(m3)), std=float(np.std(m3)))
    print("m3 mean", float(np.mean(m3)), flush=True)
    # null: identical call rates => no membership-correlated missingness
    null = [call_pattern_channel(s, missing_member=0.03, missing_non=0.03) for s in range(10)]
    res["m3_null"] = dict(per_seed=null, mean=float(np.mean(null)), std=float(np.std(null)))
    print("m3 null mean", float(np.mean(null)), flush=True)
    # matched no-holdout variant (self-reference bias of the missingness channel)
    m3_selfref = [call_pattern_channel(s, holdout=False) for s in range(5)]
    res["m3_selfref"] = dict(mean=float(np.mean(m3_selfref)), std=float(np.std(m3_selfref)))
    print("m3 self-ref mean", float(np.mean(m3_selfref)), flush=True)
    # interpretable magnitude sweep
    res["m3_sweep"] = m3_sweep()
    print("m3 sweep", json.dumps(res["m3_sweep"]), flush=True)
    # (2) structure
    print("=== structure ===", flush=True)
    sim_fst = []
    for s in range(10):
        pop = simulate_population(n_groups=3, n_per_group=450, n_snps=1363, ld_blocks=100, seed=s)
        sim_fst.append(fst_like(pop["G"], pop["ancestry"]))
    real_fst = {}
    for c in (1, 2, 22):
        try:
            G, sup, pos = get_real_panel_chr(c, max_snps=2000)  # 2000 markers, fast
            real_fst[f"chr{c}"] = fst_like(G, sup)
        except FileNotFoundError:
            pass
    res["structure"] = dict(sim_fst_measured=[float(np.mean(sim_fst)), float(np.std(sim_fst))],
                            sim_fst_nominal=0.08, real_fst=real_fst)
    print("structure", json.dumps(res["structure"]), flush=True)
    with open(os.path.join(OUT, "s20_extra.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
