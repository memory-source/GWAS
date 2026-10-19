"""Li-Stephens semi-synthetic cohort from the real 1000G phased library.

The real-data arm of the paper is capped at the 1000G sample set (~3300
haplotypes); the paper cannot tell whether the real PC<AF ordering is real
structure or small n. This module builds large cohorts that carry REAL haplotype
structure by mosaic haplotype copying (Li-Stephens / HAPGEN2-style):

  - per super-population, the donor library is the real phased haplotypes
    (H0 and H1 rows) of that super-population;
  - each synthetic individual's two haplotypes are independent mosaic copies:
    walk the markers, and at marker j switch to a fresh random donor haplotype
    with probability 1 - exp(-rho * gap_bp), where gap_bp is the physical gap
    from the previous marker;
  - the switch points therefore follow the standard HMM; the copied segments
    inherit the donor's real allele values, hence the real LD structure.

Membership labels are assigned AFTER generation by the same 3-way disjoint split
used everywhere else, so there is no label circularity.

Controls (mandatory): LD decay r^2(d) of the cohort vs the donor library;
relatedness (pairwise longest-shared-run) of the cohort vs the donor library;
AF spectrum; top-10 eigenvalue fraction. `controls_report` returns all four.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LD_EDGES = (0, 2000, 10000, 50000, 200000)


def load_donors(chrom=22, stride_to=20000):
    """Real per-super-population donor library, strided to <= stride_to markers.

    Returns dict(lib (pop -> stacked H0/H1 int8 (2n_p, M)), dosage (pop -> real
    individual dosages int8 (n_p, M)), pos, cm, M)."""
    from real_phased import get_real_phased, SUPERPOP
    import genetic_map
    H0, H1, sup, pos = get_real_phased(chrom, ld_window=1, max_snps=200000)
    m = H0.shape[1]
    keep = (np.linspace(0, m - 1, min(stride_to, m)).astype(int)
            if m > stride_to else np.arange(m))
    H0, H1, pos = H0[:, keep], H1[:, keep], pos[keep]
    cm = genetic_map.bp_to_cm(pos, chrom)
    lib, dosage = {}, {}
    for p in SUPERPOP:
        idx = np.where(sup == p)[0]
        lib[p] = np.vstack([H0[idx], H1[idx]]).astype(np.int8)
        dosage[p] = (H0[idx] + H1[idx]).astype(np.int8)
    return dict(lib=lib, dosage=dosage, pos=pos, cm=cm, M=int(len(keep)))


def _switch_prob(pos, rho):
    """Per-marker switch probability 1 - exp(-rho * gap_bp); first marker 0."""
    gap = np.diff(pos, prepend=pos[0]).astype(np.float64)
    gap[0] = 0.0
    return 1.0 - np.exp(-rho * gap)


def copy_block(donors, n, psw, rng, block=4000):
    """Generate n mosaic haplotypes from `donors` (D, M). Returns (n, M) int8.

    Vectorized per block: a candidate donor index and an independent switch coin
    are drawn per (individual, marker); the ACTIVE donor at marker j is the
    candidate at the last switch position <= j (a running maximum of the switch
    index), which realizes the Li-Stephens HMM exactly."""
    D, M = donors.shape
    donorsT = np.ascontiguousarray(donors.T)                 # (M, D)
    out = np.empty((n, M), dtype=np.int8)
    cols = np.arange(M)[None, :]
    for s in range(0, n, block):
        e = min(s + block, n)
        b = e - s
        cand = rng.integers(0, D, size=(b, M))
        sw = rng.random((b, M)) < psw[None, :]
        sw[:, 0] = True
        idx = np.where(sw, cols, -1)
        last = np.maximum.accumulate(idx, axis=1)            # last switch col <= j
        active = np.take_along_axis(cand, last, axis=1)      # (b, M) donor index
        out[s:e] = donorsT[cols, active].astype(np.int8)
    return out


def generate_cohort(chrom=22, n_per_group=1000, rho=1e-8, seed=0, stride_to=20000,
                    return_hap=False):
    """3*n_per_group individuals per super-population (enough for the 3-way
    split), mosaic-copied from the real donor library. With return_hap=True the
    ACTUAL generated mosaic haplotypes (H0/H1) are returned alongside the dosage,
    so a relatedness control can read the true copied phase rather than a dosage
    phase proxy."""
    from real_phased import SUPERPOP
    d = load_donors(chrom, stride_to=stride_to)
    lib, pos = d["lib"], d["pos"]
    psw = _switch_prob(pos, rho)
    rng = np.random.default_rng(seed)
    pool = 3 * n_per_group
    Gs, anc, pops, H0s, H1s = [], [], [], [], []
    for k, p in enumerate(SUPERPOP):
        h0 = copy_block(lib[p], pool, psw, rng)
        h1 = copy_block(lib[p], pool, psw, rng)
        Gs.append((h0 + h1).astype(np.int8))
        anc.append(np.full(pool, k))
        pops.append(np.full(pool, p, dtype=object))
        if return_hap:
            H0s.append(h0); H1s.append(h1)
    G = np.concatenate(Gs, axis=0)
    out = dict(G=G, ancestry=np.concatenate(anc).astype(np.int64),
               pop=np.concatenate(pops), n_groups=3, pos=pos, cm=d["cm"],
               rho=rho, chrom=chrom, M=d["M"])
    if return_hap:
        out["H0"] = np.concatenate(H0s, axis=0)
        out["H1"] = np.concatenate(H1s, axis=0)
    return out


# ----------------------------------------------------------------- controls ---
def ld_decay(G, pos, edges=LD_EDGES, seed=0, n_pairs=8000):
    """Mean r^2 in physical-distance bins (bp) from a dosage matrix (n, M)."""
    rng = np.random.default_rng(seed)
    m = len(pos)
    j = rng.integers(0, m, size=n_pairs)
    step = rng.integers(1, max(2, m // 4), size=n_pairs)
    k = np.minimum(j + step, m - 1)
    ok = j != k
    j, k = j[ok], k[ok]
    a = G[:, j].astype(np.float64); b = G[:, k].astype(np.float64)
    sel = (a.std(0) > 0) & (b.std(0) > 0)
    a, b = a[:, sel], b[:, sel]
    aa = a - a.mean(0); bb = b - b.mean(0)
    r = (aa * bb).sum(0) / np.sqrt((aa ** 2).sum(0) * (bb ** 2).sum(0))
    r2 = r ** 2
    dbp = np.abs(pos[k[sel]].astype(np.float64) - pos[j[sel]].astype(np.float64))
    curve, centers = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (dbp >= lo) & (dbp < hi)
        curve.append(float(r2[s].mean()) if s.sum() else float("nan"))
        centers.append(float(dbp[s].mean()) if s.sum() else float("nan"))
    return curve, centers


def relatedness(H, n_target=150, n_study=200, seed=0):
    """Mean, over target haplotypes, of the longest exact-match run (markers) to
    the best-matching OTHER haplotype within a study subset of H — the pairwise
    IBD proxy. (Full-marker per target, so it is exact; a small study subset
    keeps the cost bounded.)"""
    rng = np.random.default_rng(seed)
    h, m = H.shape
    ti = rng.permutation(h)[:n_target]
    si = rng.permutation(h)[:min(n_study, h)]
    S = H[si].astype(np.int16)
    idx = np.arange(m, dtype=np.int32)[None, :]
    pos_in_s = -np.ones(h, dtype=np.int64)
    pos_in_s[si] = np.arange(len(si))
    best = np.zeros(len(ti))
    for a, t in enumerate(ti):
        E = S == H[t].astype(np.int16)[None, :]
        sidx = np.minimum(np.maximum.accumulate(np.where(E, 0, idx + 1), axis=1), m - 1)
        span = np.where(E, idx - sidx + 1, 0)
        rb = span.max(1)
        sp = pos_in_s[t]
        if sp >= 0:
            rb[sp] = 0
        best[a] = rb.max()
    return float(best.mean()), float(best.std())


def af_spectrum(G, edges=(0.05, 0.15, 0.3, 0.5)):
    af = G.mean(0) / 2.0
    maf = np.minimum(af, 1 - af)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (maf >= lo) & (maf < hi)
        out.append(float(maf[s].mean()) if s.sum() else float("nan"))
    return out


def parsimony_phase(G):
    """Real-individual phase proxy: het sites alternate 0/1 along the marker
    order. Returns a haplotype matrix (2n, M) int8 from a dosage matrix (n, M)."""
    G = np.asarray(G)
    het = G == 1
    alt = np.cumsum(het, axis=1) % 2          # alternating parity at het sites
    h0 = np.where(het, alt, G // 2).astype(np.int8)
    h1 = (G - h0).astype(np.int8)
    return np.vstack([h0, h1])


def controls_report(cohort, donors, seed=0):
    """All four controls: LD decay, relatedness, AF spectrum, top-10 eig frac.
    Relatedness is computed on the ACTUAL copied haplotype phase when the cohort
    carries it (generate_cohort(..., return_hap=True)); otherwise it falls back
    to a dosage parsimony-phase proxy, which understates the true copied-haplotype
    sharing and is labeled as a proxy in the report."""
    from real_phased import SUPERPOP
    G = cohort["G"]; pos = cohort["pos"]
    Dd = np.concatenate([donors["dosage"][p] for p in SUPERPOP], axis=0)
    coh_ld, centers = ld_decay(G, pos, seed=seed)
    don_ld, _ = ld_decay(Dd, pos, seed=seed)
    if "H0" in cohort and "H1" in cohort:
        Htrue = np.vstack([cohort["H0"], cohort["H1"]])
        coh_rel = relatedness(Htrue, seed=seed)
        rel_phase = "actual"
        don_lib = np.concatenate([donors["lib"][p] for p in SUPERPOP], axis=0)
        don_rel = relatedness(don_lib, seed=seed)
    else:
        coh_rel = relatedness(parsimony_phase(G), seed=seed)
        don_rel = relatedness(parsimony_phase(Dd), seed=seed)
        rel_phase = "parsimony_proxy"

    def eig_frac(X):
        Xc = X.astype(np.float64); Xc = Xc - Xc.mean(0, keepdims=True)
        s = np.linalg.svd(Xc, compute_uv=False)
        ev = s ** 2
        return float(ev[:10].sum() / ev.sum())
    return dict(ld_cohort=coh_ld, ld_donor=don_ld, ld_centers_bp=centers,
                relatedness_cohort=coh_rel, relatedness_donor=don_rel,
                relatedness_phase=rel_phase,
                af_spectrum_cohort=af_spectrum(G), af_spectrum_donor=af_spectrum(Dd),
                eig10_frac_cohort=eig_frac(G), eig10_frac_donor=eig_frac(Dd),
                M=int(cohort["M"]), n_per_group=int(G.shape[0] // 3))


if __name__ == "__main__":
    import json, time
    t0 = time.time()
    d = load_donors(22, stride_to=20000)
    print("donors", {p: v.shape for p, v in d["lib"].items()}, "M", d["M"],
          "span_bp", int(d["pos"][0]), int(d["pos"][-1]), flush=True)
    for rho in (1e-9, 2e-9, 5e-9, 1e-8):
        c = generate_cohort(22, n_per_group=150, rho=rho, seed=0, stride_to=20000)
        rep = controls_report(c, d)
        print("rho", rho, "ld_coh", [round(x, 3) for x in rep["ld_cohort"]],
              "ld_don", [round(x, 3) for x in rep["ld_donor"]],
              "rel", round(rep["relatedness_cohort"][0], 1),
              round(rep["relatedness_donor"][0], 1),
              "eig10", round(rep["eig10_frac_cohort"], 3),
              round(rep["eig10_frac_donor"], 3), flush=True)
    print("elapsed", round(time.time() - t0, 1), "s", flush=True)
