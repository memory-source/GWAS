"""S32 - Real PHASED (H0/H1) IBD attack (B5).

S26 measured a dosage-level IBS window-sharing attack; the paper explicitly
defers "exploiting H0/H1 separately". This closes that gap with a genuinely
phased long-IBD attack on real 1000G:

  - Phased haplotype arrays H0/H1 (already cached by real_phased.py).
  - Attack: count, per target haplotype and per target, the maximum length of a
    run of consecutive markers that exactly matches ANY study haplotype
    (longest shared run, LSR). This is the standard IBD-segment proxy; it uses
    phase, which the dosage IBS score cannot.
  - Positive control: PLANT mathematical IBD by copying a donor haplotype
    segment into a member subset -> the attack must detect it (invariant).
  - Real arm: members are a random subset of the same super-pops (unrelated, no
    private IBD) -> LSR should be <= the AF channel, bounding the phased
    haplotype channel among unrelated individuals.
  - Length thresholds L in {20, 50, 100} markers control the segment length.

Output: results/raw/s32_phased_ibd.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from common import af_lrt_attack, auroc

OUT = "/data/lvvyongji1/GWAS/results/raw"


def _runs_maxlen(match):
    """Given a boolean (n_hap, m) match matrix, return per-row longest run of True."""
    n, m = match.shape
    best = np.zeros(n, dtype=np.int32)
    cur = np.zeros(n, dtype=np.int32)
    for j in range(m):
        cur = np.where(match[:, j], cur + 1, 0)
        best = np.maximum(best, cur)
    return best


def lsr_score(target_H0, target_H1, study_H, window=None, max_markers=None):
    """Longest exact shared run (in markers) of a target haplotype with ANY
    study haplotype, per target individual (max over H0/H1 and over study haps).

    target_H0/H1: (T, m); study_H: (S, m) with S = 2*n_study haplotypes."""
    m = study_H.shape[1]
    if max_markers is not None:
        m = min(m, max_markers)
    out = np.zeros(target_H0.shape[0], dtype=np.int32)
    # chunk over haplotypes to bound memory (matches broadcast is T x S x m)
    S = study_H.shape[0]
    for targ in (target_H0, target_H1):
        T = targ.shape[0]
        best = np.zeros(T, dtype=np.int32)
        for s0 in range(0, S, 50):
            blk = study_H[s0:s0 + 50, :m]                    # (b, m)
            eq = (targ[:, :m, None] == blk[None, :, :].transpose(0, 2, 1))  # (T, m, b)
            best = np.maximum(best, _runs_maxlen(eq.any(2)))
        out = np.maximum(out, best)
    return out


def plant_phased_ibd(H0, H1, sup, members_frac=0.3, seg=150, seed=0):
    """Copy a donor's H0 segment into a member subset (per superpop) -> planted
    haplotype IBD detectable by LSR."""
    rng = np.random.default_rng(seed)
    H0 = H0.copy()
    m = H0.shape[1]
    start = int(rng.integers(0, m - seg))
    for p in np.unique(sup):
        idx = np.where(sup == p)[0]
        donor = idx[0]
        n_mem = int(len(idx) * members_frac)
        for i in idx[1:1 + n_mem]:
            H0[i, start:start + seg] = H0[donor, start:start + seg]
    return H0, (start, seg)


def split_idx(sup, sups, per=150, seed=0):
    rng = np.random.default_rng(seed)
    si, ri, ni = [], [], []
    for p in sups:
        idx = rng.permutation(np.where(sup == p)[0])
        q = min(per, idx.size // 3)
        si.append(idx[:q]); ri.append(idx[q:2 * q]); ni.append(idx[2 * q:3 * q])
    return np.concatenate(si), np.concatenate(ri), np.concatenate(ni)


def positive_control(chrom=22, per=150, n_splits=5):
    from real_phased import get_real_phased, SUPERPOP
    H0, H1, sup, pos = get_real_phased(chrom, ld_window=1000, max_snps=20000)
    rows = []
    for s in range(n_splits):
        P0, seg = plant_phased_ibd(H0, H1, sup, seed=s)
        si, ri, ni = split_idx(sup, SUPERPOP, per=per, seed=s)
        studyH = np.vstack([P0[si], H1[si]])                 # (2*|si|, m)
        # exclude self-match: a member trivially matches its own hap, so score
        # each member (and each non-member) against study haps with that
        # individual's OWN two rows removed from the study set.
        def lsr_excl(inds, memb):
            sc = np.zeros(len(inds), dtype=np.int32)
            for t, ix in enumerate(inds):
                keep = np.ones(studyH.shape[0], bool)
                keep[t] = False; keep[len(inds) + t] = False
                sc[t] = lsr_score(P0[ix:ix + 1], H1[ix:ix + 1], studyH[keep],
                                  max_markers=4000)[0]
            return sc
        sm_ex = lsr_excl(si, True)
        sn_ex = lsr_excl(ni, False)
        y = np.concatenate([np.ones(len(si)), np.zeros(len(ni))])
        planted = float(auroc(np.concatenate([sm_ex, sn_ex]), y))
        rows.append(dict(split=s, planted_lsr=planted, seg=seg,
                         mean_mem=float(sm_ex.mean()), mean_non=float(sn_ex.mean())))
    agg = {k: [float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows]))]
           for k in ("planted_lsr",)}
    return dict(agg=agg, per_split=rows, n_snps=int(H0.shape[1]))


def real_arm(chrom=22, ld_window=1500, per=150, n_splits=8, L=None):
    from real_phased import get_real_phased, SUPERPOP
    L = L or [20, 50, 100]
    H0, H1, sup, pos = get_real_phased(chrom, ld_window=ld_window, max_snps=40000)
    res = {}
    for thr in L:
        rows = []
        for s in range(n_splits):
            si, ri, ni = split_idx(sup, SUPERPOP, per=per, seed=s)
            studyH = np.vstack([H0[si], H1[si]])
            sm = lsr_score(H0[si], H1[si], studyH)
            sn = lsr_score(H0[ni], H1[ni], studyH)
            # AUC of "LSR >= thr" as a classifier (segment presence)
            y = np.concatenate([np.ones(len(si)), np.zeros(len(ni))])
            score = np.concatenate([sm, sn]).astype(float)
            lsr_auc = float(auroc(score, y))
            hit = (np.concatenate([sm, sn]) >= thr).astype(float)
            hit_auc = float(auroc(hit, y))
            # AF baseline
            siteG = H0[si] + H1[si]; refG = H0[ri] + H1[ri]; nonG = H0[ni] + H1[ni]
            af_s, af_r = siteG.mean(0) / 2.0, refG.mean(0) / 2.0
            am = af_lrt_attack(siteG, af_s, af_r); an = af_lrt_attack(nonG, af_s, af_r)
            af = float(auroc(np.concatenate([am, an]), y))
            rows.append(dict(split=s, lsr_auc=lsr_auc, hit_auc=hit_auc, af=af,
                             mean_lsr_mem=float(sm.mean()), mean_lsr_non=float(sn.mean())))
        res[f"L{thr}"] = {k: [float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows]))]
                          for k in ("lsr_auc", "hit_auc", "af", "mean_lsr_mem", "mean_lsr_non")}
        print("real chr%d L=%d" % (chrom, thr), json.dumps(
            {k: round(v[0], 4) for k, v in res[f"L{thr}"].items()}), flush=True)
    return dict(thresholds=L, per_L=res, n_snps=int(H0.shape[1]),
                n_study_haps=int(2 * per))


def run():
    res = {}
    t0 = time.time()
    try:
        res["positive_control"] = positive_control()
        print("planted LSR", res["positive_control"]["agg"], flush=True)
    except Exception as e:
        res["positive_control"] = dict(error=str(e))
    try:
        res["real"] = real_arm()
    except Exception as e:
        res["real"] = dict(error=str(e))
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "s32_phased_ibd.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("done", res["elapsed_s"], "s", flush=True)


if __name__ == "__main__":
    run()
