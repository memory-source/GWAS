"""Phased 1000G loader: keeps haplotype phase for a genuine IBD/haplotype attack.

Fact-check: the 1000G Phase-3 VCFs in data/external are PHASED (GT fields use
'|', e.g. `0|1`; filenames are `...shapeit2_mvncall_integrated...`). The existing
parsers (real_panel.py, real_multi.py) collapse GT to dosage and discard phase.
This module keeps the two haplotypes separately, which enables a genuine
phased-haplotype / long-identity-by-state (IBS) membership attack.

Cache: results/cache/real_phased_chr<C>_<tag>.npz with H0, H1 (n x m int8),
sup, pos.
"""
import gzip
import hashlib
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXT = os.path.join(ROOT, "data", "external")
CACHE = os.path.join(ROOT, "results", "cache")
PANEL = os.path.join(EXT, "1000g_panel.txt")
SUPERPOP = ("AFR", "EAS", "EUR")


def _vcf_path(chrom):
    for c in (f"chr{chrom}.vcf.gz",
              f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
              f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"):
        p = os.path.join(EXT, c)
        if os.path.exists(p) and os.path.getsize(p) > 10_000:
            return p
    raise FileNotFoundError(f"no VCF for chr{chrom}")


def load_panel():
    samples, sup = [], []
    with open(PANEL) as f:
        next(f)
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                samples.append(parts[0]); sup.append(parts[2])
    return samples, np.array(sup)


def _parse_phased(chrom, max_snps=80000, maf_lo=0.05):
    vcf = _vcf_path(chrom)
    samples, sup = load_panel()
    keep = np.array([s in SUPERPOP for s in sup])
    ids = [s for s, k in zip(samples, keep) if k]
    sample_sup = sup[keep]
    H0_rows, H1_rows, pos, cols = [], [], [], None
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
                c = {h: i for i, h in enumerate(hdr)}
                cols = [c[s] for s in ids]
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(cols):
                continue
            ref, alt = p[3], p[4]
            if "," in alt or "." in alt or len(ref) != 1 or len(alt) != 1:
                continue
            h0 = np.empty(len(cols), dtype=np.int8)
            h1 = np.empty(len(cols), dtype=np.int8)
            bad = False
            for j, i in enumerate(cols):
                g = p[i]
                a = g[0]
                sep = g[1] if len(g) > 1 else "/"
                b = g[2] if len(g) > 2 else "0"
                if a == "." or b == ".":
                    bad = True; break
                h0[j] = 1 if a == "1" else 0
                h1[j] = 1 if b == "1" else 0
            if bad:
                continue
            af = (h0.mean() + h1.mean()) / 2.0
            maf = min(af, 1 - af)
            if maf < maf_lo:
                continue
            H0_rows.append(h0); H1_rows.append(h1); pos.append(int(p[1]))
            if len(H0_rows) >= max_snps:
                break
            if len(H0_rows) % 20000 == 0:
                print(f"  chr{chrom} phased: {len(H0_rows)}", flush=True)
    H0 = np.array(H0_rows, dtype=np.int8).T
    H1 = np.array(H1_rows, dtype=np.int8).T
    return H0, H1, sample_sup, np.array(pos)


def get_real_phased(chrom, ld_window=1000, max_snps=80000, use_cache=True):
    os.makedirs(CACHE, exist_ok=True)
    tag = hashlib.md5(f"{chrom}_{max_snps}_{ld_window}".encode()).hexdigest()[:8]
    path = os.path.join(CACHE, f"real_phased_chr{chrom}_{tag}.npz")
    if use_cache and os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        return z["H0"], z["H1"], z["sup"], z["pos"]
    H0, H1, sup, pos = _parse_phased(chrom, max_snps=max_snps)
    keep, last = [], -10 ** 12
    for j in range(len(pos)):
        if pos[j] - last >= ld_window:
            keep.append(j); last = pos[j]
    keep = np.array(keep)
    H0, H1, pos = H0[:, keep], H1[:, keep], pos[keep]
    np.savez(path, H0=H0, H1=H1, sup=sup, pos=pos)
    return H0, H1, sup, pos


if __name__ == "__main__":
    for c in (22,):
        H0, H1, sup, pos = get_real_phased(c)
        print(f"chr{c} phased", H0.shape, H1.shape, {p: int((sup == p).sum()) for p in SUPERPOP})
