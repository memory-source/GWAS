"""S41 direct self-match: the released haplotype set trivially leaks membership.

The s41 LSR statistic removes each candidate's OWN two haplotypes before scoring
(a leave-self-out relatedness diagnostic). But the study RELEASES its haplotype
set, which for every member literally contains that member's two haplotypes: a
candidate who holds the release and its own genotypes finds itself verbatim. This
probe measures that trivial membership leak directly — score = number of the
candidate's own haplotypes present verbatim in the released set — and reports it
next to the leave-self-out LSR statistic, so the two are not conflated.

Output: results/raw/s41_directmatch.json
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, "/data/lvvyongji1/GWAS/code")
from common import auroc
from real_phased import get_real_phased, SUPERPOP

OUT = "/data/lvvyongji1/GWAS/results/raw"


def one(chrom=22, per=80, seed=0, max_snps=50000):
    H0, H1, sup, pos = get_real_phased(chrom, ld_window=1, max_snps=200000)
    m = H0.shape[1]
    keep = np.linspace(0, m - 1, max_snps).astype(int) if m > max_snps else np.arange(m)
    H0, H1 = H0[:, keep], H1[:, keep]
    rng = np.random.default_rng(seed)
    si, ri, ni = [], [], []
    for p in SUPERPOP:
        idx = rng.permutation(np.where(sup == p)[0])
        q = min(per, idx.size // 3)
        si.append(idx[:q]); ri.append(idx[q:2 * q]); ni.append(idx[2 * q:3 * q])
    si, ni = np.concatenate(si), np.concatenate(ni)
    released = np.vstack([H0[si], H1[si]])           # the released haplotype set
    keys = {released[i].tobytes() for i in range(released.shape[0])}

    def self_match(h0, h1):
        return np.array([int(h0[i].tobytes() in keys) + int(h1[i].tobytes() in keys)
                         for i in range(len(h0))], float)

    sm = self_match(H0[si], H1[si])
    sn = self_match(H0[ni], H1[ni])
    y = np.concatenate([np.ones(len(si)), np.zeros(len(ni))])
    return dict(n_members=int(len(si)), n_non=int(len(ni)),
                self_match_auc=float(auroc(np.concatenate([sm, sn]), y)),
                frac_members_present=float((sm > 0).mean()),
                frac_non_present=float((sn > 0).mean()))


def run():
    out = {}
    for chrom in (22, 1, 2):
        out[f"chr{chrom}"] = one(chrom)
        print(f"chr{chrom}", json.dumps({k: (round(v, 4) if isinstance(v, float) else v)
              for k, v in out[f"chr{chrom}"].items()}), flush=True)
    with open(os.path.join(OUT, "s41_directmatch.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out


if __name__ == "__main__":
    run()
