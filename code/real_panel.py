"""Shared cache for the parsed 1000 Genomes chr22 panel.

The VCF stream-parse in s7_real is the expensive step; every downstream
real-data experiment (orientation study, matched-operating-point arm) reuses
the same G/pos/sup arrays. Cache to results/cache/real_panel.npz keyed by
the parse parameters so a parameter change invalidates cleanly.
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
VCF = os.path.join(EXT, "chr22.vcf.gz")
PANEL = os.path.join(EXT, "1000g_panel.txt")
SUPERPOP = ("AFR", "EAS", "EUR")


def load_panel():
    samples, sup = [], []
    with open(PANEL) as f:
        next(f)
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            samples.append(parts[0])
            sup.append(parts[2])
    return samples, np.array(sup)


def _parse_vcf(max_snps=60000, maf_lo=0.05, maf_hi=0.95):
    samples, sup = load_panel()
    keep = np.array([s in SUPERPOP for s in sup])
    sample_ids = [s for s, k in zip(samples, keep) if k]
    sample_sup = sup[keep]
    G_rows, pos, sample_cols = [], [], None
    with gzip.open(VCF, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
                cols = {h: i for i, h in enumerate(hdr)}
                sample_cols = [cols[s] for s in sample_ids]
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(sample_cols):
                continue
            ref, alt = p[3], p[4]
            if "," in alt or "." in alt or len(ref) != 1 or len(alt) != 1:
                continue
            row = np.empty(len(sample_cols), dtype=np.int8)
            bad = False
            for j, i in enumerate(sample_cols):
                g = p[i]
                a = g[0]
                b = g[2] if len(g) > 2 else "0"
                if a == "." or b == ".":
                    bad = True
                    break
                row[j] = (1 if a == "1" else 0) + (1 if b == "1" else 0)
            if bad:
                continue
            af = row.mean() / 2.0
            maf = min(af, 1 - af)
            if maf < maf_lo or maf > maf_hi:
                continue
            G_rows.append(row)
            pos.append(int(p[1]))
            if len(G_rows) >= max_snps:
                break
    G = np.array(G_rows, dtype=np.int8)
    assert G.size and G.max() > 1
    return G.T, sample_sup, np.array(pos)


def get_real_panel(ld_window=15000, max_snps=60000, use_cache=True):
    """Return (G, sup, pos) MAF-filtered + distance-thinned, cached to npz."""
    os.makedirs(CACHE, exist_ok=True)
    tag = hashlib.md5(f"{max_snps}_{ld_window}".encode()).hexdigest()[:8]
    path = os.path.join(CACHE, f"real_panel_{tag}.npz")
    if use_cache and os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        return z["G"], z["sup"], z["pos"]
    G, sup, pos = _parse_vcf(max_snps=max_snps)
    keep = []
    last = -10 ** 12
    for j in range(len(pos)):
        if pos[j] - last >= ld_window:
            keep.append(j)
            last = pos[j]
    keep = np.array(keep)
    G, pos = G[:, keep], pos[keep]
    np.savez(path, G=G, sup=sup, pos=pos)
    return G, sup, pos


def three_way_split(sup, per=150, seed=0):
    """Disjoint study/reference/non-member per super-population."""
    rng = np.random.default_rng(seed)
    sidx, ridx, nidx = [], [], []
    for p in SUPERPOP:
        idx = rng.permutation(np.where(sup == p)[0])
        per_p = min(per, idx.size // 3)
        sidx.append(idx[:per_p])
        ridx.append(idx[per_p:2 * per_p])
        nidx.append(idx[2 * per_p:3 * per_p])
    return np.concatenate(sidx), np.concatenate(ridx), np.concatenate(nidx)


if __name__ == "__main__":
    G, sup, pos = get_real_panel()
    print("real panel", G.shape, "pos", pos[0], pos[-1],
          {p: int((sup == p).sum()) for p in SUPERPOP})
