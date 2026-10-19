"""S36 - Phased IBD with self-exclusion (fixes the s32 self-match leak).

S32's real-arm LSR saturates at the panel length for members because a member's
own haplotype is IN the study set -> it matches itself end-to-end. The correct
attack excludes the target individual's own two haplotypes from the study set
(an attacker scoring a candidate does not count the candidate as a reference).
This block re-runs BOTH the planted positive control and the real arm with that
exclusion, on chr22 (and chr1/chr2 if their phased caches exist).

Output: results/raw/s36_phased_ibd.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from common import af_lrt_attack, auroc
from s32_phased_ibd import lsr_score, plant_phased_ibd, split_idx

OUT = "/data/lvvyongji1/GWAS/results/raw"


def lsr_selfexcl(target_idx, H0, H1, study_idx, studyH, m=None):
    """Per-target LSR against study haps EXCLUDING that target's own two rows.

    Requires target_idx to be a subset of study_idx (same ordering prefix used to
    build studyH), so the self rows are at positions t and |study_idx|+t."""
    T = len(target_idx)
    n_st = len(study_idx)
    sc = np.zeros(T, dtype=np.int32)
    for t, ix in enumerate(target_idx):
        keep = np.ones(studyH.shape[0], bool)
        keep[t] = False; keep[n_st + t] = False
        sc[t] = lsr_score(H0[ix:ix + 1], H1[ix:ix + 1], studyH[keep], max_markers=m)[0]
    return sc


def positive_control(chrom=22, per=150, n_splits=4):
    from real_phased import get_real_phased, SUPERPOP
    H0, H1, sup, pos = get_real_phased(chrom, ld_window=1000, max_snps=20000)
    rows = []
    for s in range(n_splits):
        P0, seg = plant_phased_ibd(H0, H1, sup, seed=s)
        si, ri, ni = split_idx(sup, SUPERPOP, per=per, seed=s)
        blockH = np.vstack([P0[si], H1[si]])
        sm = lsr_selfexcl(si, P0, H1, si, blockH)
        # non-members are NOT in the study set -> no self rows to remove
        sn = lsr_score(P0[ni], H1[ni], blockH)
        y = np.concatenate([np.ones(len(si)), np.zeros(len(ni))])
        planted = float(auroc(np.concatenate([sm, sn]), y))
        # unplanted control: same split on the ORIGINAL haplotypes
        blockH0 = np.vstack([H0[si], H1[si]])
        sm0 = lsr_selfexcl(si, H0, H1, si, blockH0)
        sn0 = lsr_score(H0[ni], H1[ni], blockH0)
        unplanted = float(auroc(np.concatenate([sm0, sn0]), y))
        rows.append(dict(split=s, planted=planted, unplanted=unplanted, seg=seg))
        print("split%d planted=%.3f unplanted=%.3f" % (s, planted, unplanted), flush=True)
    agg = {k: [float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows]))]
           for k in ("planted", "unplanted")}
    return dict(agg=agg, per_split=rows, n_snps=int(H0.shape[1]))


def real_arm(chrom=22, per=150, n_splits=6):
    from real_phased import get_real_phased, SUPERPOP
    H0, H1, sup, pos = get_real_phased(chrom, ld_window=1500, max_snps=40000)
    rows = []
    for s in range(n_splits):
        si, ri, ni = split_idx(sup, SUPERPOP, per=per, seed=s)
        studyH = np.vstack([H0[si], H1[si]])
        sm = lsr_selfexcl(si, H0, H1, si, studyH)
        sn = lsr_score(H0[ni], H1[ni], studyH)
        score = np.concatenate([sm, sn]).astype(float)
        y = np.concatenate([np.ones(len(si)), np.zeros(len(ni))])
        lsr_auc = float(auroc(score, y))
        # AF baseline
        siteG = H0[si] + H1[si]; refG = H0[ri] + H1[ri]; nonG = H0[ni] + H1[ni]
        am = af_lrt_attack(siteG, siteG.mean(0) / 2.0, refG.mean(0) / 2.0)
        an = af_lrt_attack(nonG, siteG.mean(0) / 2.0, refG.mean(0) / 2.0)
        af = float(auroc(np.concatenate([am, an]), y))
        rows.append(dict(split=s, lsr_auc=lsr_auc, af=af,
                         mean_mem=float(sm.mean()), mean_non=float(sn.mean()),
                         panel=int(H0.shape[1])))
        print("real split%d lsr=%.3f af=%.3f mem=%.0f non=%.0f" % (
            s, lsr_auc, af, sm.mean(), sn.mean()), flush=True)
    agg = {k: [float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows]))]
           for k in ("lsr_auc", "af", "mean_mem", "mean_non")}
    return dict(agg=agg, per_split=rows, n_snps=int(H0.shape[1]))


def run():
    res = {}
    t0 = time.time()
    try:
        res["positive_control"] = positive_control()
        print("planted agg", res["positive_control"]["agg"], flush=True)
    except Exception as e:
        import traceback; traceback.print_exc()
        res["positive_control"] = dict(error=str(e))
    try:
        res["real_chr22"] = real_arm()
    except Exception as e:
        import traceback; traceback.print_exc()
        res["real_chr22"] = dict(error=str(e))
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "s36_phased_ibd.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("done", res["elapsed_s"], "s", flush=True)


if __name__ == "__main__":
    run()
