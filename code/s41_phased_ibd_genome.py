"""S41 - multi-chromosome physical-length phased-IBD bound.

Extends s36 (chr22 only, run length in MARKERS, single density) to chr1+chr2+
chr22 with the IBD statistic measured in physical length (centiMorgans, via the
1000G Phase-3 IMPUTE2 genetic map). Two arms per chromosome, identical marker set:

  - PLANTED control: a donor haplotype segment of ~10 cM is copied into a member
    subset; the statistic must detect those members (calibrates the scale).
  - REAL arm: members are a random subset of the same super-populations
    (unrelated, no private IBD); each target's own two released haplotypes are
    removed before scoring, so a member can only be detected through a REAL long
    match to OTHER released individuals. The bound tested: real IBD <= AF.

Statistic = per target haplotype, the longest run of consecutive markers that
matches a SINGLE released haplotype exactly (max over released haplotypes),
converted to its centiMorgan span. This is the pairwise IBD-segment proxy from
s36, now physical-length. Length-threshold sweep {0.5,1,2,3,5} cM: the fraction
of targets carrying a segment >= threshold, and the presence-indicator AUROC.

Output: results/raw/s41_phased_ibd_genome.json
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from common import auroc, af_lrt_attack
import genetic_map

OUT = "/data/lvvyongji1/GWAS/results/raw"
CHROMS = (1, 2, 22)
MAX_MARKERS = 50000            # per chromosome, uniform stride (documented below)
THRESHOLDS_CM = (0.5, 1.0, 2.0, 3.0, 5.0)
PLANT_CM = 10.0


def load_chrom(chrom, stride_to=MAX_MARKERS):
    """Full-density phased cache, uniformly strided to <= stride_to markers so
    every chromosome gets genome-wide coverage at a reported density, plus the
    cumulative cM at each retained marker."""
    from real_phased import get_real_phased
    H0, H1, sup, pos = get_real_phased(chrom, ld_window=1, max_snps=200000)
    m = H0.shape[1]
    if m > stride_to:
        keep = np.linspace(0, m - 1, stride_to).astype(int)
    else:
        keep = np.arange(m)
    H0, H1, pos = H0[:, keep], H1[:, keep], pos[keep]
    cm = genetic_map.bp_to_cm(pos, chrom)
    return dict(H0=H0, H1=H1, sup=sup, pos=pos, cm=cm, n_markers=int(len(keep)))


def _cm_run_vs_one(targets_H, ref_hap, cm):
    """Per target haplotype, the longest exact-match run to `ref_hap` measured in
    cM (the cM span is cm[end]-cm[start] for the run ending at each matching
    cell; its row maximum is the longest cM segment, which need not be the run
    with the most markers when marker spacing varies)."""
    E = targets_H == ref_hap[None, :]
    m = targets_H.shape[1]
    idx = np.arange(m, dtype=np.int32)[None, :]
    sidx = np.maximum.accumulate(np.where(E, 0, idx + 1), axis=1)
    # clip the run-start index: when the last marker mismatches, `sidx` reaches m.
    # cm[start] is only read on matching cells, so a clip to m-1 is harmless there
    # and removes the out-of-bounds read on the cells where E is False.
    span = np.where(E, cm[None, :] - cm[np.minimum(sidx, m - 1)], -np.inf)
    return span.max(axis=1)


def best_cm_per_target(targets_H, released_H, cm, n_members):
    """Per target individual, the best exact-match cM run against released_H,
    excluding the target's OWN two released haplotypes.

    targets_H: (2T, m) = [hap0 of every target; hap1 of every target], targets
    ordered members-then-non-members (T = n_members + n_non). released_H:
    (2*n_members, m) = [members' hap0; members' hap1]. released row j<n_members
    belongs to member j's hap0 (must be excluded from target rows j), and
    j>=n_members to member j-n_members's hap1 (excluded from target row
    T+(j-n_members)). A target's own reference haplotypes are skipped for THAT
    comparison (masked before the running max), so earlier legitimate matches
    are preserved, and the returned per-individual score is the max over the
    individual's two (hap0/hap1) target rows."""
    T2 = targets_H.shape[0]
    T = T2 // 2
    accum = np.full(T2, -np.inf)
    for j in range(released_H.shape[0]):
        cmrun = _cm_run_vs_one(targets_H, released_H[j], cm)
        if j < n_members:
            cmrun[j] = -np.inf                       # member j's own hap0 row
            cmrun[T + j] = -np.inf                   # member j's own hap1 row
        else:
            mj = j - n_members
            cmrun[mj] = -np.inf                      # member mj's own hap0 row
            cmrun[T + mj] = -np.inf                  # member mj's own hap1 row
        accum = np.maximum(accum, cmrun)
    return accum


def split_idx(sup, sups, per, seed):
    rng = np.random.default_rng(seed)
    si, ri, ni = [], [], []
    for p in sups:
        idx = rng.permutation(np.where(sup == p)[0])
        q = min(per, idx.size // 3)
        si.append(idx[:q]); ri.append(idx[q:2 * q]); ni.append(idx[2 * q:3 * q])
    return np.concatenate(si), np.concatenate(ri), np.concatenate(ni)


def plant_segment(H0, sup, si, coord, seg_cM=PLANT_CM, members_frac=0.5, seed=0):
    """Copy a donor's H0 segment spanning ~seg_cM into a subset of the STUDY
    MEMBERS (`si`), per super-population. Only members receive the segment, so
    the statistic must separate them from non-members. Returns the modified H0 and
    the realized cM span."""
    rng = np.random.default_rng(seed)
    H0 = H0.copy()
    cm = coord["cm"]
    total = cm[-1] - cm[0]
    start_cm = float(rng.uniform(cm[0] + 0.1 * total, cm[-1] - 0.1 * total - seg_cM))
    lo = int(np.searchsorted(cm, start_cm))
    hi = max(int(np.searchsorted(cm, start_cm + seg_cM)), lo + 5)
    for p in np.unique(sup):
        pm = si[sup[si] == p]
        if len(pm) < 2:
            continue
        donor = int(pm[0])
        n_mem = max(1, int(len(pm) * members_frac))
        recips = pm[1:1 + n_mem]
        H0[np.ix_(recips, np.arange(lo, hi))] = H0[donor, lo:hi]
    return H0, float(cm[min(hi, len(cm) - 1)] - cm[lo])


def af_baseline(H0, H1, si, ri, ni, mask=None):
    siteG = H0[si] + H1[si]; refG = H0[ri] + H1[ri]; nonG = H0[ni] + H1[ni]
    af_s = siteG.mean(0) / 2.0; af_r = refG.mean(0) / 2.0
    am = af_lrt_attack(siteG, af_s, af_r, mask=mask)
    an = af_lrt_attack(nonG, af_s, af_r, mask=mask)
    y = np.concatenate([np.ones(len(si)), np.zeros(len(ni))])
    return float(auroc(np.concatenate([am, an]), y))


def arm(coord, si, ri, ni, planted, seed):
    """One chromosome, one split: real and planted LSR-cM vs AF.

    Target layout: rows [0:T) = hap0 of every target, rows [T:2T) = hap1, with
    targets ordered members-then-non-members (T = n_members + n_non)."""
    H0 = coord["H0"]; H1 = coord["H1"]; cm = coord["cm"]; sup = coord["sup"]
    n_st = len(si); n_non = len(ni); T = n_st + n_non
    hap0 = np.vstack([H0[si], H0[ni]])          # (T, m) hap0 of members then non
    hap1 = np.vstack([H1[si], H1[ni]])          # (T, m) hap1
    TH = np.vstack([hap0, hap1])                # (2T, m)
    released = np.vstack([H0[si], H1[si]])
    score_hap = best_cm_per_target(TH, released, cm, n_members=n_st)
    # per-INDIVIDUAL score = max over the individual's two haplotypes
    score_cm = np.maximum(score_hap[:T], score_hap[T:])          # (T,)
    mem_cm = score_cm[:n_st]; non_cm = score_cm[n_st:T]
    y = np.concatenate([np.ones(n_st), np.zeros(n_non)])
    out = dict(lsr_auc=float(auroc(np.concatenate([mem_cm, non_cm]), y)),
               mean_cm_mem=float(mem_cm.mean()), mean_cm_non=float(non_cm.mean()),
               max_cm_mem=float(mem_cm.max()))
    for thr in THRESHOLDS_CM:
        hit = np.concatenate([mem_cm >= thr, non_cm >= thr]).astype(float)
        out[f"auc_cm{thr}"] = float(auroc(hit, y))
        out[f"frac_mem_cm{thr}"] = float((mem_cm >= thr).mean())
        out[f"frac_non_cm{thr}"] = float((non_cm >= thr).mean())
    out["af"] = af_baseline(H0, H1, si, ri, ni)
    if planted:
        P0, seg_cM = plant_segment(coord["H0"], sup, si, coord, seed=seed)
        # planted member-recipients carry a copied donor segment; non-members do
        # not. Score against the planted released set with self-exclusion.
        relP = np.vstack([P0[si], H1[si]])
        hap0p = np.vstack([P0[si], P0[ni]]); hap1p = np.vstack([H1[si], H1[ni]])
        THp = np.vstack([hap0p, hap1p])
        scp = best_cm_per_target(THp, relP, cm, n_members=n_st)
        scp = np.maximum(scp[:T], scp[T:])
        mp = scp[:n_st]; npn = scp[n_st:T]
        out["planted_lsr_auc"] = float(auroc(np.concatenate([mp, npn]), y))
        out["planted_seg_cM"] = seg_cM
        out["planted_mean_cm_mem"] = float(mp.mean())
        out["planted_mean_cm_non"] = float(npn.mean())
        for thr in THRESHOLDS_CM:
            hitp = np.concatenate([mp >= thr, npn >= thr]).astype(float)
            out[f"planted_auc_cm{thr}"] = float(auroc(hitp, y))
    return out


def run(chrom=22, per=80, n_splits=3, do_planted=True, stride_to=MAX_MARKERS):
    from real_phased import SUPERPOP
    coord = load_chrom(chrom, stride_to=stride_to)
    rows = []
    for s in range(n_splits):
        si, ri, ni = split_idx(coord["sup"], SUPERPOP, per=per, seed=s)
        r = arm(coord, si, ri, ni, planted=do_planted, seed=s)
        rows.append(r)
        print("chr%d split%d lsr=%.3f planted=%.3f af=%.3f mean_cm(m/n)=%.2f/%.2f" % (
            chrom, s, r["lsr_auc"], r.get("planted_lsr_auc", float("nan")), r["af"],
            r["mean_cm_mem"], r["mean_cm_non"]), flush=True)
    keys = [k for k in rows[0] if isinstance(rows[0][k], float)]
    agg = {k: [float(np.mean([r[k] for r in rows])), float(np.std([r[k] for r in rows]))]
           for k in keys}
    return dict(chrom=chrom, agg=agg, per_split=rows,
                n_markers=coord["n_markers"], span_bp=[int(coord["pos"][0]), int(coord["pos"][-1])],
                span_cM=float(coord["cm"][-1] - coord["cm"][0]))


def run_all():
    res = {"per_chrom": {}, "thresholds_cM": list(THRESHOLDS_CM), "per": 80,
           "n_splits": 3, "plant_cM": PLANT_CM}
    t0 = time.time()
    for c in CHROMS:
        try:
            res["per_chrom"][f"chr{c}"] = run(chrom=c)
        except FileNotFoundError as e:
            res["per_chrom"][f"chr{c}"] = dict(error=str(e))
            print("skip chr", c, e, flush=True)
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "s41_phased_ibd_genome.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("done", res["elapsed_s"], "s", flush=True)
    return res


if __name__ == "__main__":
    run_all()
