"""S26 - genuine haplotype-sharing / IBS membership attack.

The paper's haplotype channel (s19) is an adjacent-PAIR joint-table score fit
in-sample on the study, so its 0.99 is a fitted-table artifact (it transfers at
chance). This adds a distinct segment-based attack:

  - an 8-marker WINDOW-SHARING count: the fraction of windows a target matches
    to any STUDY individual vs any REFERENCE individual (unphased dosages).

Real arm runs on unphased dosages; real_phased.py loads the 1000G haplotypes
but this arm does not yet exploit H0/H1 separately.

Controls: a planted-IBD positive control (copy a segment into a member subset ->
the attack MUST detect it) and the real-1000G arm.

Output: results/raw/s26_ibd.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import simulate_population, split_population3, allele_freq, auroc

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "raw")


def window_match_score(target_G, study_G, ref_G, window=8, max_mismatch=0,
                       exclude_self=False):
    """Vectorized haplotype-sharing score over non-overlapping windows.

    For each target and each window of `window` markers, ask whether the
    target's window matches ANY study row within `max_mismatch` mismatches, and
    the same against the reference panel. Score = (#study-matched windows) -
    (#ref-matched windows). An exact IBS/haplotype-sharing proxy, vectorized
    over windows (no per-marker Python loop)."""
    m = study_G.shape[1]
    nw = m // window
    T = target_G.shape[0]
    tgt = target_G[:, :nw * window].reshape(T, nw, window)
    stu = study_G[:, :nw * window].reshape(study_G.shape[0], nw, window)
    rf = ref_G[:, :nw * window].reshape(ref_G.shape[0], nw, window)
    sm = np.zeros(T); sr = np.zeros(T)
    self_idx = np.arange(T)
    for w in range(nw):
        tw = tgt[:, w, :]                                   # (T, window)
        ds = (tw[:, None, :] != stu[None, :, w, :]).sum(-1)  # (T, S)
        if exclude_self:
            ds[self_idx, self_idx] = 10 ** 6
        sm += (ds <= max_mismatch).any(1)
        dr = (tw[:, None, :] != rf[None, :, w, :]).sum(-1)   # (T, R)
        sr += (dr <= max_mismatch).any(1)
    return sm - sr


def ibs_score(target_G, study_G, ref_G, max_mismatch=0, cap=None,
              exclude_self_from_study=False, window=8):
    return window_match_score(target_G, study_G, ref_G, window=window,
                              max_mismatch=max_mismatch,
                              exclude_self=exclude_self_from_study)


def plant_ibd(pop, members_frac=0.3, seg=200, seed=0):
    """Copy a contiguous segment from a donor (index 0) into a subset of members,
    creating a private IBD segment the attack should detect."""
    rng = np.random.default_rng(seed)
    G = pop["G"].copy()
    anc = pop["ancestry"]
    m = G.shape[1]
    start = rng.integers(0, m - seg)
    for k in np.unique(anc):
        idx = np.where(anc == k)[0]
        n_mem = int(len(idx) * members_frac)
        donor = idx[0]
        for i in idx[1:1 + n_mem]:
            G[i, start:start + seg] = G[donor, start:start + seg]
    pop = dict(pop); pop["G"] = G; pop["ibd_segment"] = (int(start), int(seg))
    return pop


def sim_positive_control(seeds=(0, 1, 2, 3, 4)):
    rows = []
    for seed in seeds:
        pop = simulate_population(n_groups=3, n_per_group=1000, n_snps=4000, seed=seed)
        pop_ibd = plant_ibd(pop, seed=seed)
        site, ref, non = split_population3(pop_ibd, 300, 300, 300)
        sm = ibs_score(site["G"], site["G"], ref["G"], max_mismatch=0,
                       exclude_self_from_study=True)
        sn = ibs_score(non["G"], site["G"], ref["G"], max_mismatch=0)
        y = np.concatenate([np.ones(len(sm)), np.zeros(len(sn))])
        planted = float(auroc(np.concatenate([sm, sn]), y))
        # negative: same but no planted IBD
        site0, ref0, non0 = split_population3(pop, 300, 300, 300)
        sm0 = ibs_score(site0["G"], site0["G"], ref0["G"], max_mismatch=0,
                        exclude_self_from_study=True)
        sn0 = ibs_score(non0["G"], site0["G"], ref0["G"], max_mismatch=0)
        unplanted = float(auroc(np.concatenate([sm0, sn0]), y))
        rows.append(dict(seed=seed, planted_ibs=planted, unplanted_ibs=unplanted))
        print("sim", json.dumps(rows[-1]), flush=True)
    agg = {k: [float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows]))]
           for k in ("planted_ibs", "unplanted_ibs")}
    return dict(per_seed=rows, agg=agg)


def real_phased_arm(chroms=(22,), per=150, n_splits=10):
    """Phased-haplotype IBS attack on real 1000G. Members are a random subset of
    the same super-populations, so there is no private IBD -> the attack should
    be near chance / <= the real AF channel."""
    from real_phased import get_real_phased, SUPERPOP
    res = {}
    for c in chroms:
        try:
            H0, H1, sup, pos = get_real_phased(c, ld_window=1000, max_snps=20000)
        except FileNotFoundError as e:
            print("skip chr", c, e, flush=True); continue
        rows = []
        for s in range(n_splits):
            rng = np.random.default_rng(s)
            sidx, ridx, nidx = [], [], []
            for p in SUPERPOP:
                idx = rng.permutation(np.where(sup == p)[0])
                q = min(per, idx.size // 3)
                sidx.append(idx[:q]); ridx.append(idx[q:2 * q]); nidx.append(idx[2 * q:3 * q])
            si, ri, ni = np.concatenate(sidx), np.concatenate(ridx), np.concatenate(nidx)
            # collapse to dosage for the phasing-free proxy, and to phased for the real one
            siteG = (H0[si] + H1[si]); refG = (H0[ri] + H1[ri]); nonG = (H0[ni] + H1[ni])
            y = np.concatenate([np.ones(len(si)), np.zeros(len(ni))])
            sm = ibs_score(siteG, siteG, refG, exclude_self_from_study=True)
            sn = ibs_score(nonG, siteG, refG)
            ibs = float(auroc(np.concatenate([sm, sn]), y))
            # AF baseline on the same split
            af_s, af_r = siteG.mean(0) / 2.0, refG.mean(0) / 2.0
            from common import af_lrt_attack
            am = af_lrt_attack(siteG, af_s, af_r); an = af_lrt_attack(nonG, af_s, af_r)
            af = float(auroc(np.concatenate([am, an]), y))
            rows.append(dict(split=s, ibs=ibs, af=af))
        agg = {k: [float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows]))]
               for k in ("ibs", "af")}
        res[f"chr{c}"] = dict(agg=agg, per_split=rows, n_snps=int(H0.shape[1]))
        print("real", c, json.dumps({k: round(v[0], 4) for k, v in agg.items()}), flush=True)
    return res


def run():
    res = dict(sim_control=sim_positive_control(), real=real_phased_arm())
    with open(os.path.join(OUT, "s26_ibd.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
