"""S36 - Phased release channel test (corrected).

Question: does releasing PHASED haplotypes (H0/H1) create a membership signal
beyond the dosage AF channel?

A prior design scored members' longest-shared-run against a study set that
INCLUDED their own haplotypes (self-match; a design flaw, not a leak). This
version instead:

  * PLANTED control: copy a donor's H0 segment into a member subset. The attack
    must detect those members (near-perfect self-match to the donor's copied
    segment) -> validates the statistic.
  * REAL arm: a candidate's own two haplotypes are removed from the released
    set before scoring, so a member can only be detected via a *real* long
    match to OTHER released individuals.

Design constraint discovered while building this: a run of >= L markers matching
ANY of 2*n_study released haplotypes is not a rare event even for unrelated
individuals, because the max over many haplotypes accumulates matches. So the
correct long-range statistic is the max over released haplotypes of the longest
run matching a SINGLE released haplotype (pairwise), i.e. the classic IBD proxy,
not "any". We therefore compute, for each target haplotype, the longest run of
consecutive markers that matches a SINGLE released haplotype exactly (max over
released haplotypes), which is only large for genuine IBD / near-duplicates.

HMM: pairwise exact match over 30 markers at ~65% per-marker match prob has
probability ~0.65^30 ~ 2e-6 per window, so unrelated individuals rarely reach
L=30; a planted copy reaches the full segment length (~300).

Output: results/raw/s36_phased_ibd.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from common import auroc

OUT = "/data/lvvyongji1/GWAS/results/raw"


def _max_run_per_row(match):
    """match: boolean (n, m) -> per-row longest run of True (no Python-per-marker)."""
    n, m = match.shape
    out = np.zeros(n, dtype=np.int32)
    for i in range(n):
        row = match[i]
        if not row.any():
            continue
        d = np.diff(np.concatenate([[0], row.astype(np.int8), [0]]))
        st = np.where(d == 1)[0]; en = np.where(d == -1)[0]
        out[i] = int((en - st).max())
    return out


def pairwise_longest_run(target_H, released_H):
    """For each target haplotype, the longest run matching a SINGLE released hap
    (max over released). target_H (T, m), released_H (S, m). -> (T,)."""
    T, m = target_H.shape
    best = np.zeros(T, dtype=np.int32)
    for s in range(released_H.shape[0]):
        eq = target_H == released_H[s][None, :]            # (T, m) exact match to hap s
        best = np.maximum(best, _max_run_per_row(eq))
    return best


def topk_meantop(target_H, released_H, k=5):
    """For each target hap, mean of the k longest runs over a SINGLE released hap.
    A smoother version of pairwise_longest_run."""
    T, m = target_H.shape
    store = np.zeros((T, k), dtype=np.int32)
    for s in range(released_H.shape[0]):
        r = _max_run_per_row(target_H == released_H[s][None, :])
        # insert r into the running top-k (vectorized-ish)
        combined = np.concatenate([store, r[:, None]], axis=1)
        store = np.sort(combined, axis=1)[:, -k:]
    return store.mean(1)


def score_targets(t_idx, H0, H1, released_H, self_blocks, method="topk", k=5):
    """Score each target individual (its two haps) against released_H, excluding
    that individual's own two haplotypes. self_blocks: array (n_target, 2) of
    released-row indices that are the target's own haps."""
    T = len(t_idx)
    out = np.zeros(T)
    for t in range(T):
        sb = self_blocks[t]
        if sb is None:
            R = released_H                       # non-member: no self rows to drop
        else:
            keep = np.ones(released_H.shape[0], bool)
            keep[sb] = False
            R = released_H[keep]
        a = pairwise_longest_run(H0[t_idx[t]:t_idx[t] + 1], R)[0]
        b = pairwise_longest_run(H1[t_idx[t]:t_idx[t] + 1], R)[0]
        out[t] = max(a, b)
    return out.astype(float)


def run(chrom=22, per=80, n_splits=5, seg=300):
    from real_phased import get_real_phased, SUPERPOP
    from s32_phased_ibd import plant_phased_ibd
    H0, H1, sup, pos = get_real_phased(chrom, ld_window=1500, max_snps=30000)
    out = {"real": [], "planted": []}
    for s in range(n_splits):
        rng = np.random.default_rng(s)
        si, ni = [], []
        for p in SUPERPOP:
            idx = rng.permutation(np.where(sup == p)[0])
            q = min(per, idx.size // 3)
            si.append(idx[:q]); ni.append(idx[2 * q:3 * q])
        si, ni = np.concatenate(si), np.concatenate(ni)
        n_st = len(si)
        released = np.vstack([H0[si], H1[si]])               # (2 n_st, m)
        self_mem = [np.array([t, n_st + t]) for t in range(n_st)]   # own two rows
        # non-members are NOT in the released set -> no self rows to remove (None)
        self_non = [None] * len(ni)
        # ---- REAL arm
        sm = score_targets(si, H0, H1, released, self_mem)
        sn = score_targets(ni, H0, H1, released, self_non)
        y = np.concatenate([np.ones(len(si)), np.zeros(len(ni))])
        auc_real = float(auroc(np.concatenate([sm, sn]), y))
        out["real"].append(dict(split=s, auc=auc_real, mem=float(sm.mean()), non=float(sn.mean())))
        # ---- PLANTED control (donor segment copied into a member subset)
        P0, pseg = plant_phased_ibd(H0, H1, sup, members_frac=0.6, seg=seg, seed=s)
        relP = np.vstack([P0[si], H1[si]])
        smP = score_targets(si, P0, H1, relP, self_mem)
        snP = score_targets(ni, H0, H1, relP, self_non)
        auc_plant = float(auroc(np.concatenate([smP, snP]), y))
        out["planted"].append(dict(split=s, auc=auc_plant, seg=pseg,
                                   mem=float(smP.mean()), non=float(snP.mean())))
        print("split%d real=%.3f planted=%.3f seg=%s" % (s, auc_real, auc_plant, pseg), flush=True)
    agg = {k: {kk: [float(np.mean([r[kk] for r in v])), float(np.std([r[kk] for r in v]))]
               for kk in ("auc", "mem", "non")} for k, v in out.items()}
    return dict(agg=agg, per_split=out, n_snps=int(H0.shape[1]), per=per)


if __name__ == "__main__":
    t0 = time.time()
    res = run()
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "s36_phased_ibd.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("agg planted", res["agg"]["planted"], flush=True)
    print("agg real", res["agg"]["real"], flush=True)
    print("done", res["elapsed_s"], "s", flush=True)
