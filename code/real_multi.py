"""Multi-chromosome 1000 Genomes loader (extends real_panel to chr1/chr2/chr22).

The chr22-only version is kept for backward compatibility; this module parses an
arbitrary 1000G Phase-3 autosomal chromosome with the same MAF filter +
distance-thinning, so the real-data studies can be replicated on independent
chromosomes (chr1, chr2) and aggregated.

Cache: results/cache/real_panel_chr<C>_<tag>.npz
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
SUPERPOP5 = ("AFR", "AMR", "EAS", "EUR", "SAS")


def _vcf_path(chrom):
    """Resolve the VCF filename for a chromosome across the known 1000G naming."""
    cands = [
        f"chr{chrom}.vcf.gz",
        f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    ]
    for c in cands:
        p = os.path.join(EXT, c)
        if os.path.exists(p) and os.path.getsize(p) > 10_000:
            return p
    raise FileNotFoundError(f"no VCF for chr{chrom} in {EXT} (tried {cands})")


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


def _parse_vcf(chrom, max_snps=200000, maf_lo=0.05, maf_hi=0.95, sups=None):
    vcf = _vcf_path(chrom)
    samples, sup = load_panel()
    sups = tuple(sups) if sups else SUPERPOP
    keep = np.array([s in sups for s in sup])
    sample_ids = [s for s, k in zip(samples, keep) if k]
    sample_sup = sup[keep]
    G_rows, pos, sample_cols = [], [], None
    with gzip.open(vcf, "rt") as f:
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
            if len(G_rows) % 20000 == 0:
                print(f"  chr{chrom}: parsed {len(G_rows)}", flush=True)
    G = np.array(G_rows, dtype=np.int8)
    assert G.size and G.max() > 1, f"chr{chrom} GT parsing observed no variation"
    return G.T, sample_sup, np.array(pos)


def get_real_panel_chr(chrom, ld_window=15000, max_snps=200000, use_cache=True, sups=None):
    os.makedirs(CACHE, exist_ok=True)
    tag = hashlib.md5(f"{chrom}_{max_snps}_{ld_window}_{sups}".encode()).hexdigest()[:8]
    path = os.path.join(CACHE, f"real_panel_chr{chrom}_{tag}.npz")
    if use_cache and os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        return z["G"], z["sup"], z["pos"]
    G, sup, pos = _parse_vcf(chrom, max_snps=max_snps, sups=sups)
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


def get_multi_chrom(chroms=(1, 2), ld_window=15000, max_snps=200000, sups=None):
    """Concatenate MD-filtered + thinned panels across chromosomes (markers
    stacked; samples aligned because the 1000G release uses a fixed sample set)."""
    Gs, poss, sup = [], [], None
    for c in chroms:
        G, s, p = get_real_panel_chr(c, ld_window=ld_window, max_snps=max_snps, sups=sups)
        Gs.append(G); poss.append(p); sup = s
    G = np.concatenate(Gs, axis=1)
    pos = np.concatenate(poss)
    return G, sup, pos


def three_way_split(sup, per=150, seed=0):
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
    for c in (1, 2):
        try:
            G, sup, pos = get_real_panel_chr(c)
            print(f"chr{c} panel", G.shape, "pos", pos[0], pos[-1],
                  {p: int((sup == p).sum()) for p in SUPERPOP})
        except FileNotFoundError as e:
            print("skip", c, e)
