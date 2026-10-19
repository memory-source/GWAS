"""S52 - second, independent haplotype channel estimator (round-15 closure).

The paper's haplotype-channel negative is currently only for the longest-shared-
run / IBS-window statistic (s36). A reviewer can object that the negative is one
estimator. This adds two competing haplotype/IBD statistics on the same real 1000G
phased data, both leave-self-out (a member's own two haplotypes are removed from
the released set before it is scored), plus the same planted positive control:

  - `max_ibs_share`: the max over released haplotypes of the whole-haplotype
    per-marker identical-by-state (IBS) share of a target haplotype (a smooth
    whole-haplotype similarity, not a run length);
  - `n_close_haps` : the number of released haplotypes whose IBS share to the
    target exceeds a fixed threshold (a count statistic robust to a single
    spurious match).

Success: both are near chance on the real unrelated arm and at least one detects
the planted segment (the planted control copies one donor's H0 segment of length
>=20% of the marker panel into the member subset, enough to clear the
max-over-all-released-haplotypes noise floor); then the haplotype negative is not
an artifact of one estimator.

Output: results/raw/s52_haplo_second_chr<C>.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from common import auroc

OUT = "/data/lvvyongji1/GWAS/results/raw"
THR = 0.78


def max_ibs_share(target_H, released_H):
    """(T,) max over released haps of the per-marker IBS share."""
    best = np.zeros(target_H.shape[0])
    for s in range(released_H.shape[0]):
        share = (target_H == released_H[s][None, :]).mean(1)
        best = np.maximum(best, share)
    return best


def n_close_haps(target_H, released_H, thr=THR):
    """(T,) count of released haps with IBS share > thr."""
    cnt = np.zeros(target_H.shape[0])
    for s in range(released_H.shape[0]):
        share = (target_H == released_H[s][None, :]).mean(1)
        cnt += (share > thr)
    return cnt


def score_targets(t_idx, H0, H1, released_H, self_blocks):
    T = len(t_idx)
    out = np.zeros((T, 2))                      # (max_ibs_share, n_close), best of 2 haps
    for t in range(T):
        sb = self_blocks[t]
        R = released_H if sb is None else released_H[np.setdiff1d(np.arange(released_H.shape[0]), sb)]
        a = max(0.0, max_ibs_share(H0[t_idx[t]:t_idx[t] + 1], R)[0],
                max_ibs_share(H1[t_idx[t]:t_idx[t] + 1], R)[0])
        c = max(n_close_haps(H0[t_idx[t]:t_idx[t] + 1], R)[0],
                n_close_haps(H1[t_idx[t]:t_idx[t] + 1], R)[0])
        out[t] = (a, c)
    return out


def run(chrom=22, per=80, n_splits=5, seg=1500):
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
        rel = np.vstack([H0[si], H1[si]])
        self_mem = [np.array([t, n_st + t]) for t in range(n_st)]
        self_non = [None] * len(ni)
        y = np.concatenate([np.ones(len(si)), np.zeros(len(ni))])
        sm = score_targets(si, H0, H1, rel, self_mem)
        sn = score_targets(ni, H0, H1, rel, self_non)
        au_real = [float(auroc(np.concatenate([sm[:, c], sn[:, c]]), y)) for c in range(2)]
        # PLANTED control: copy a donor H0 segment into the STUDY-MEMBER subset
        # only (si), at a length >=20% of the marker panel so the IBS-share lift
        # clears the max-over-all-released-haps noise floor.
        P0_si, pseg = plant_phased_ibd(H0[si], H1[si], sup[si], members_frac=0.6,
                                       seg=seg, seed=s)
        relP = np.vstack([P0_si, H1[si]])
        smP = score_targets(np.arange(n_st), P0_si, H1[si], relP, self_mem)
        snP = score_targets(np.arange(len(ni)), H0[ni], H1[ni], relP, self_non)
        au_plant = [float(auroc(np.concatenate([smP[:, c], snP[:, c]]), y)) for c in range(2)]
        out["real"].append(dict(split=s, auc_max_ibs=au_real[0], auc_n_close=au_real[1]))
        out["planted"].append(dict(split=s, auc_max_ibs=au_plant[0], auc_n_close=au_plant[1]))
        print("split%d real maxIBS=%.3f nclose=%.3f | planted %.3f/%.3f seg=%d" % (
            s, au_real[0], au_real[1], au_plant[0], au_plant[1], pseg[1]), flush=True)
    agg = {}
    for k, v in out.items():
        agg[k] = {}
        for kk in ("auc_max_ibs", "auc_n_close"):
            vals = [r[kk] for r in v]
            agg[k][kk] = [float(np.mean(vals)), float(np.std(vals))]
    return dict(agg=agg, per_split=out, n_snps=int(H0.shape[1]), per=per, thr=THR, seg=seg)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22)
    a = ap.parse_args()
    t0 = time.time()
    res = run(chrom=a.chrom)
    res["chrom"] = a.chrom
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"s52_haplo_second_chr{a.chrom}.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("agg real", res["agg"]["real"], flush=True)
    print("agg planted", res["agg"]["planted"], flush=True)
    print("done", res["elapsed_s"], "s", flush=True)
