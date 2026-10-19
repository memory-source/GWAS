"""S7 - Real-data external validity on 1000 Genomes chr22.

Parses the public Phase-3 chr22 VCF (GT field), selects three super-population
pseudo-sites (EUR/EAS/AFR), splits each three ways DISJOINTLY (study members /
attacker reference / independent non-members), MAF-filters and LD-thins, then
reruns the frontier + PC-spectrum mechanism on real released allele frequencies.

Output: results/raw/s7_real.json
"""
import gzip
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import allele_freq, af_lrt_attack, anc_only_attack, auroc, tpr_at_fpr, pca

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXT = os.path.join(ROOT, "data", "external")
OUT = os.path.join(ROOT, "results", "raw")
VCF = os.path.join(EXT, "chr22.vcf.gz")
PANEL = os.path.join(EXT, "1000g_panel.txt")

SUPERPOP = {"EUR", "EAS", "AFR"}


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


def parse_vcf(max_snps=60000, maf_lo=0.05, maf_hi=0.95):
    """Stream the VCF once; keep biallelic SNVs passing MAF, up to max_snps.

    NOTE: GT is parsed from the fully split line (the sample block is indexed by
    absolute column position), and BOTH alleles are read, so genotypes are in
    {0,1,2}. A previous version split only 9 fields, leaving p[9] as the entire
    sample block as one string, and never read the second allele.
    """
    samples, sup = load_panel()
    keep_samples = np.array([s in SUPERPOP for s in sup])
    sample_ids = [s for s, k in zip(samples, keep_samples) if k]
    sample_sup = sup[keep_samples]
    print("kept samples:", len(sample_ids),
          {p: int((sample_sup == p).sum()) for p in SUPERPOP}, flush=True)

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
            chrom, pos_, _id, ref, alt = p[0], p[1], p[2], p[3], p[4]
            if "," in alt or "." in alt or len(ref) != 1 or len(alt) != 1:
                continue  # biallelic SNV only
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
            pos.append(int(pos_))
            if len(G_rows) >= max_snps:
                break
            if len(G_rows) % 5000 == 0:
                print("parsed", len(G_rows), flush=True)
    G = np.array(G_rows, dtype=np.int8)          # (n_snps, n_samples)
    print("final G", G.shape, flush=True)
    assert G.size and G.max() > 1, "GT parsing did not observe any het/hom-alt"
    return G.T, sample_sup, np.array(pos)          # (n_samples, n_snps)


def ld_thin(G, pos, max_snps=8000, window=15_000):
    """Distance-thinning (>=window bp apart), a cheap linkage proxy, NOT r^2
    pruning. Named accordingly."""
    keep = []
    last = -10**12
    for j in range(len(pos)):
        if pos[j] - last >= window:
            keep.append(j)
            last = pos[j]
        if len(keep) >= max_snps:
            break
    return np.array(keep)


def three_way_split(G, sup, per=150, seed=0):
    """Disjoint study/reference/non-member panels per super-population.

    Within each super-population the sample order is SHUFFLED with a recorded
    seed before splitting, because the 1000G panel file groups sub-populations;
    contiguous thirds would confound the split with sub-population composition.
    """
    rng = np.random.default_rng(seed)
    sidx, ridx, nidx = [], [], []
    for p in sorted(SUPERPOP):
        idx = rng.permutation(np.where(sup == p)[0])
        per_p = min(per, idx.size // 3)
        sidx.append(idx[:per_p])
        ridx.append(idx[per_p:2 * per_p])
        nidx.append(idx[2 * per_p:3 * per_p])
    return np.concatenate(sidx), np.concatenate(ridx), np.concatenate(nidx)


def boot(scores, labels, fpr=0.01, n=300, seed=0):
    rng = np.random.default_rng(seed)
    pos = np.where(labels == 1)[0]; neg = np.where(labels == 0)[0]
    v = [tpr_at_fpr(scores[np.r_[rng.choice(pos, pos.size, True), rng.choice(neg, neg.size, True)]],
                    labels[np.r_[rng.choice(pos, pos.size, True), rng.choice(neg, neg.size, True)]], fpr)
         for _ in range(n)]
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def run(split_seeds=(0, 1, 2, 3, 4)):
    """Average the frontier across several random 3-way splits of the real
    panels (split seed only reshuffles assignment; the genotypes are fixed)."""
    G, sup, pos = parse_vcf()
    keep = ld_thin(G, pos)
    G = G[:, keep]
    print("after LD-thin", G.shape, flush=True)

    tier_auroc = {}
    tier_tpr = {}
    ctrl_struct, ctrl_pc = [], []
    for s in split_seeds:
        sidx, ridx, nidx = three_way_split(G, sup, seed=s)
        site = G[sidx]; ref = G[ridx]; non = G[nidx]
        print("split", s, site.shape, flush=True)
        af_s = allele_freq(site); af_r = allele_freq(ref)
        for tier in ("F2_topk", "F4_round1", "F3_round2", "F5_full", "F5_perpop"):
            if tier == "F2_topk":
                score = np.abs(af_s - af_r)
                vis = np.zeros(G.shape[1], bool); vis[np.argsort(score)[-200:]] = True
                sm = af_lrt_attack(site, af_s, af_r, mask=vis)
                sn = af_lrt_attack(non, af_s, af_r, mask=vis)
            elif tier == "F4_round1":
                sm = af_lrt_attack(site, np.round(af_s, 1), af_r)
                sn = af_lrt_attack(non, np.round(af_s, 1), af_r)
            elif tier == "F3_round2":
                sm = af_lrt_attack(site, np.round(af_s, 2), af_r)
                sn = af_lrt_attack(non, np.round(af_s, 2), af_r)
            elif tier == "F5_perpop":
                sm = np.zeros(site.shape[0]); sn = np.zeros(non.shape[0])
                for p in sorted(SUPERPOP):
                    ms = sup[sidx] == p; mn = sup[nidx] == p; mr = sup[ridx] == p
                    afp = allele_freq(site[ms]); afrp = allele_freq(ref[mr])
                    sm[ms] = af_lrt_attack(site[ms], afp, afrp)
                    sn[mn] = af_lrt_attack(non[mn], afp, afrp)
            else:
                sm = af_lrt_attack(site, af_s, af_r)
                sn = af_lrt_attack(non, af_s, af_r)
            lab = np.concatenate([np.ones(len(sm)), np.zeros(len(sn))])
            sc = np.concatenate([sm, sn])
            tier_auroc.setdefault(tier, []).append(float(auroc(sc, lab)))
            tier_tpr.setdefault(tier, []).append(float(tpr_at_fpr(sc, lab, 0.01)))
        # controls
        _, _, load_r = pca(ref, k=10); mean_r = ref.mean(0)
        am = anc_only_attack(site, load_r, mean_r, np.ones(G.shape[1], bool))
        an = anc_only_attack(non, load_r, mean_r, np.ones(G.shape[1], bool))
        lab = np.concatenate([np.ones(len(am)), np.zeros(len(an))])
        ctrl_struct.append(float(auroc(np.concatenate([am, an]), lab)))
        _, _, load_s = pca(site, k=10); mean_s = site.mean(0)
        pm = anc_only_attack(site, load_s, mean_s, np.ones(G.shape[1], bool))
        pn = anc_only_attack(non, load_s, mean_s, np.ones(G.shape[1], bool))
        ctrl_pc.append(float(auroc(np.concatenate([pm, pn]), lab)))

    res = {}
    for tier in tier_auroc:
        res[tier] = dict(auroc=float(np.mean(tier_auroc[tier])),
                         auroc_std=float(np.std(tier_auroc[tier])),
                         tpr_at_fpr1=float(np.mean(tier_tpr[tier])),
                         tpr_std=float(np.std(tier_tpr[tier])))
        print(tier, json.dumps(res[tier]), flush=True)
    res["structure_only_control"] = float(np.mean(ctrl_struct))
    res["study_pc_distance_auroc"] = float(np.mean(ctrl_pc))
    print("controls", res["structure_only_control"], res["study_pc_distance_auroc"], flush=True)

    res["meta"] = dict(n_snps=int(G.shape[1]), n_per_panel=450,
                       split_seeds=list(split_seeds),
                       superpops=sorted(SUPERPOP))
    with open(os.path.join(OUT, "s7_real.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    run()
