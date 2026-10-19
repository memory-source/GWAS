"""Generate paper figures from results/raw/*.json into paper/figures/ (v2, uses
the disjoint-negative results and the boundary curve)."""
import glob
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "results", "raw")
FIG = os.path.join(ROOT, "paper", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.bbox": "tight"})


def load(name):
    with open(os.path.join(RAW, name)) as f:
        return json.load(f)


def fig1_frontier():
    r = load("s1c_frontier.json")["results"]
    tiers = ["F2_topk", "F4_round1", "F3_round2", "F5_full", "F5_sitewise"]
    labels = ["F2\ntop-200", "F4\nAF 1dp", "F3\nAF 2dp", "F5\nAF full", "F5\nper-site"]
    auroc = [r[t]["auroc"] for t in tiers]
    tpr = [r[t]["tpr_at_fpr1"] for t in tiers]
    lo = [r[t]["tpr_ci"][0] for t in tiers]
    hi = [r[t]["tpr_ci"][1] for t in tiers]
    anc = r["F5_full"]["ancestry_only_auroc"]
    x = np.arange(len(tiers))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].bar(x, auroc, color="#3b6ea5")
    ax[0].axhline(0.5, ls="--", c="grey", label="chance")
    ax[0].axhline(anc, ls=":", c="red", label="ancestry-only control")
    ax[0].set_xticks(x); ax[0].set_xticklabels(labels, fontsize=9)
    ax[0].set_ylabel("membership-inference AUROC"); ax[0].set_ylim(0.4, 1.0)
    ax[0].set_title("(a) Release-precision frontier"); ax[0].legend(fontsize=8)
    yerr = [np.array(tpr) - np.array(lo), np.array(hi) - np.array(tpr)]
    ax[1].bar(x, tpr, color="#a5523b", yerr=yerr, capsize=3)
    ax[1].set_xticks(x); ax[1].set_xticklabels(labels, fontsize=9)
    ax[1].set_ylabel("TPR @ FPR ≤ 1% (95% CI)")
    ax[1].set_title("(b) Low-FPR operating point")
    fig.savefig(os.path.join(FIG, "fig1_frontier.pdf"))
    fig.savefig(os.path.join(FIG, "fig1_frontier.png"))
    plt.close(fig)


def fig2_marginal():
    r = load("s2b_intermediate.json")["results"]
    keys = ["base_final_only", "base+M1_lambdaGC", "base+M4_conditionbit",
            "base+M3_missingness", "base+M2_pcaSpectrum"]
    labs = ["final only", "+λ_GC", "+cond. bit", "+missingness", "+PCA spectrum"]
    auc = [r[k]["auc_mean"] for k in keys]
    err = [r[k]["auc_std"] for k in keys]
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(np.arange(len(keys)), auc, yerr=err, color="#3b6ea5", capsize=3)
    ax.set_xticks(np.arange(len(keys))); ax.set_xticklabels(labs, rotation=15, fontsize=9)
    ax.set_ylabel("meta-attack AUROC"); ax.set_ylim(0.70, 0.90)
    ax.axhline(r["base_final_only"]["auc_mean"], ls="--", c="grey", lw=1)
    ax.annotate("+0.098 AUROC", (4, r["base+M2_pcaSpectrum"]["auc_mean"] + 0.005),
                ha="center", fontsize=9, color="#a5523b")
    ax.set_title("Marginal leakage of intermediate reveals (disjoint negatives)")
    fig.savefig(os.path.join(FIG, "fig2_marginal.pdf"))
    fig.savefig(os.path.join(FIG, "fig2_marginal.png"))
    plt.close(fig)


def fig3_boundary():
    r = load("s5_boundary.json")
    ns = sorted({v["n_site"] for v in r.values()})
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for dec, col, lab in ((1, "#a5523b", "AF rounded to 1 dp"),
                          (2, "#3b6ea5", "AF rounded to 2 dp"),
                          (None, "black", "AF full precision")):
        ys = [r[f"n{n}_dec{dec}"]["tpr_at_fpr1"] for n in ns]
        ax.plot(ns, ys, "o-", color=col, label=lab)
    ax.set_xlabel("individuals per site"); ax.set_ylabel("TPR @ FPR ≤ 1%")
    ax.set_title("Site-size / precision leakage boundary")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.savefig(os.path.join(FIG, "fig3_boundary.pdf"))
    fig.savefig(os.path.join(FIG, "fig3_boundary.png"))
    plt.close(fig)


def fig4_paradox():
    r = load("s3b_paradox.json")["summary"]
    fams = ["A_scree", "B_gap"]
    flabs = ["family A\n(scree threshold)", "family B\n(largest gap)"]
    w = [r[f + "|san=False"]["n_probes_mean"] for f in fams]
    wo = [r[f + "|san=True"]["n_probes_mean"] for f in fams]
    x = np.arange(len(fams)); bw = 0.35
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.bar(x - bw / 2, w, bw, label="no sanitized diagnostics", color="#a5523b")
    ax.bar(x + bw / 2, wo, bw, label="sanitized pre-provided", color="#3b6ea5")
    ax.set_xticks(x); ax.set_xticklabels(flabs, fontsize=9)
    ax.set_ylabel("probe reveals issued"); ax.legend(fontsize=8)
    ax.set_title("Sanitized-diagnostics paradox (scripted analyst model)")
    for i, (a, b) in enumerate(zip(w, wo)):
        if a > 0:
            ax.annotate(f"-{100*(a-b)/a:.0f}%", (i, max(a, b) + 0.1), ha="center", fontsize=9)
    fig.savefig(os.path.join(FIG, "fig4_paradox.pdf"))
    fig.savefig(os.path.join(FIG, "fig4_paradox.png"))
    plt.close(fig)


def fig5_pc_channel():
    """The PC channel: orientation, released granularity (sim vs real), and the
    PC-vs-AF crossover map. Sources: s9_orientation, s11_robust, s13_real_loadings,
    s12_crossover."""
    o = load("s9_orientation.json")
    g = load("s11_robust.json")
    gr = load("s13_real_loadings.json")
    cx = load("s12_crossover.json")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

    # (a) orientation: AF vs naive-PC vs calibrated-PC vs combined
    grp = ["AF only", "PC channel", "AF + PC"]
    sim = [o["sim_matched"]["ld_blocks=100"]["agg"][k] for k in
           ("af_only", "pc_naive", "af_plus_calibrated_pc")]
    real = [o["real"]["agg"][k] for k in
            ("af_only", "pc_naive", "af_plus_calibrated_pc")]
    x = np.arange(len(grp)); bw = 0.38
    ax[0].bar(x - bw / 2, sim, bw, label="sim (matched op. pt.)", color="#3b6ea5")
    ax[0].bar(x + bw / 2, real, bw, label="real 1000G chr22", color="#a5523b")
    ax[0].axhline(0.5, ls="--", c="grey", lw=1)
    ax[0].set_xticks(x); ax[0].set_xticklabels(grp, fontsize=8, rotation=10)
    ax[0].set_ylabel("AUROC (PC orientation calibrated)"); ax[0].set_ylim(0.0, 1.0)
    ax[0].set_title("(a) PC channel vs AF release"); ax[0].legend(fontsize=8)

    # (b) released granularity, sim vs real
    ks = ["top1_loadings_calibrated", "top2_loadings_calibrated",
          "top3_loadings_calibrated", "top5_loadings_calibrated",
          "full_loadings_calibrated", "eigenvalues_only_calibrated"]
    labs = ["top-1", "top-2", "top-3", "top-5", "full", "eig only"]
    sv = [g["s9c"]["agg"][k] for k in ks]
    rv = [gr["agg"][k] for k in ks]
    x = np.arange(len(ks))
    ax[1].plot(x, sv, "o-", color="#3b6ea5", label="sim")
    ax[1].plot(x, rv, "s-", color="#a5523b", label="real")
    ax[1].axhline(0.5, ls="--", c="grey", lw=1)
    ax[1].set_xticks(x); ax[1].set_xticklabels(labs, fontsize=9)
    ax[1].set_ylabel("AUROC (sign-calibrated)"); ax[1].set_ylim(0.4, 1.0)
    ax[1].set_title("(b) Released-PC granularity"); ax[1].legend(fontsize=8)

    # (c) PC-vs-AF increment heatmap
    nps = cx["n_per_group"]; ms = cx["n_snps"]
    M = np.full((len(nps), len(ms)), np.nan)
    for i, npg in enumerate(nps):
        for j, m in enumerate(ms):
            M[i, j] = cx["grid"][f"n{npg}_m{m}"]["delta_pc_minus_af"]
    im = ax[2].imshow(M, cmap="RdBu_r", vmin=-0.15, vmax=0.15, aspect="auto")
    for i in range(len(nps)):
        for j in range(len(ms)):
            ax[2].text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=7)
    ax[2].set_xticks(range(len(ms))); ax[2].set_xticklabels(ms)
    ax[2].set_yticks(range(len(nps)))
    ax[2].set_yticklabels([f"{n}\n({3*int(0.3*n)})" for n in nps])
    ax[2].set_xlabel("SNPs / marker count")
    ax[2].set_ylabel("population-group $n$\n(study ind.)")
    ax[2].set_title("(c) PC − AF increment (AUROC)")
    fig.colorbar(im, ax=ax[2], fraction=0.046)
    fig.savefig(os.path.join(FIG, "fig5_pc_channel.pdf"))
    fig.savefig(os.path.join(FIG, "fig5_pc_channel.png"))
    plt.close(fig)


def fig6_replication():
    """Replication panel: (a) multi-chromosome real PC-vs-AF; (b) msprime PC
    cross-check; (c) haplotype channel; (d) simulation scale sweep.
    Missing files are skipped gracefully."""
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.2))
    have = []

    def L(n):
        p = os.path.join(RAW, n)
        return load(n) if os.path.exists(p) else None

    # (a) multichrom real orientation
    r = L("s16_multichrom.json")
    if r:
        chs = list(r["per_chrom"].keys())
        af = [r["per_chrom"][c]["orientation"]["af_only"][0] for c in chs]
        pc = [r["per_chrom"][c]["orientation"]["pc_calibrated"][0] for c in chs]
        comb = [r["per_chrom"][c]["orientation"]["af_plus_pc"][0] for c in chs]
        x = np.arange(len(chs)); bw = 0.26
        ax[0].bar(x - bw, af, bw, label="AF", color="#3b6ea5")
        ax[0].bar(x, pc, bw, label="PC", color="#a5523b")
        ax[0].bar(x + bw, comb, bw, label="AF+PC", color="#5a8c5a")
        ax[0].axhline(0.5, ls="--", c="grey", lw=1)
        ax[0].set_xticks(x); ax[0].set_xticklabels(chs, fontsize=9)
        ax[0].set_ylabel("AUROC"); ax[0].set_title("(a) Real replication by chr")
        ax[0].legend(fontsize=8)
        have.append("a")

    # (b) msprime PC cross-check
    r = L("s18_msprime_pc.json")
    if r:
        a = r["agg"]
        grp = ["AF", "PC", "top1", "top3", "top5", "full", "eig only"]
        keys = ["af_only", "pc_calibrated", "top1", "top3", "top5", "full", "eig_only"]
        vals = [a[k][0] for k in keys]
        ax[1].bar(range(len(grp)), vals, color="#3b6ea5")
        ax[1].axhline(0.5, ls="--", c="grey", lw=1)
        ax[1].set_xticks(range(len(grp))); ax[1].set_xticklabels(grp, fontsize=8, rotation=15)
        ax[1].set_ylabel("AUROC"); ax[1].set_title("(b) msprime PC cross-check")
        have.append("b")

    # (c) haplotype/LD channel
    r = L("s19_haplotype.json")
    if r:
        a = r["agg"]
        grp = ["AF", "hap/LD", "hap OOS"]
        vals = [a["af_only"][0], a["hap_ld"][0], a["hap_oos_transfer"][0]]
        ax[2].bar(range(len(grp)), vals, color="#a5523b")
        ax[2].axhline(0.5, ls="--", c="grey", lw=1)
        ax[2].set_xticks(range(len(grp)))
        ax[2].set_xticklabels(grp, fontsize=8, rotation=15)
        ax[2].set_ylabel("AUROC"); ax[2].set_ylim(0.0, 1.05)
        ax[2].set_title("(c) Haplotype/LD channel")
        have.append("c")

    # (d) simulation scale sweep
    r = L("s17_scale.json")
    if r:
        cells = r["cells"]
        ns = sorted({int(k.split("_")[0][1:]) for k in cells})
        key = lambda n: f"n{n}_m10000_ld100" if f"n{n}_m10000_ld100" in cells else f"n{n}_m20000_ld100"
        af = [cells[key(n)]["af_only"] for n in ns]
        pc = [cells[key(n)]["pc_calibrated"] for n in ns]
        d = [cells[key(n)]["delta_pc_minus_af"] for n in ns]
        ax[3].plot(ns, af, "o-", color="#3b6ea5", label="AF")
        ax[3].plot(ns, pc, "s-", color="#a5523b", label="PC")
        ax[3].set_xscale("log"); ax[3].set_xlabel("study individuals")
        ax[3].set_ylabel("AUROC"); ax[3].set_title("(d) Sim scale sweep")
        ax[3].legend(fontsize=8)
        for xx, dd in zip(ns, d):
            ax[3].annotate(f"{dd:+.2f}", (xx, 0.51), ha="center", fontsize=7)
        have.append("d")

    fig.savefig(os.path.join(FIG, "fig6_replication.pdf"))
    fig.savefig(os.path.join(FIG, "fig6_replication.png"))
    plt.close(fig)


def fig7_governance():
    """Limitation-closure panel: (a) Fst crossover (ordering is one parameter);
    (b) governed vs ungoverned release; (c) cross-run composition; (d) IBD
    positive control vs real. Missing files skipped gracefully."""
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.2))

    def L(n):
        p = os.path.join(RAW, n)
        return load(n) if os.path.exists(p) else None

    # (a) Fst crossover
    r = L("s21_sensitivity.json")
    if r:
        fsts = r["fsts"]
        for npg, col in ((300, "#9bb7d4"), (800, "#3b6ea5"), (1600, "#1b3d6b")):
            ds = [r["grid"][f"f{f}_n{npg}"]["agg"]["delta_pc_minus_af"] for f in fsts]
            ax[0].plot(fsts, ds, "o-", color=col, label=f"n={npg}")
        ax[0].axhline(0, ls="--", c="grey", lw=1)
        # the crossover is size-specific. On the nominal axis the gap crosses
        # zero between Fst 0.04-0.08 at n=1600, 0.08-0.12 at n=800, and
        # 0.16-0.24 at n=300; mark one band per site size.
        for lo, hi, col in ((0.04, 0.08, "#1b3d6b"),
                            (0.08, 0.12, "#3b6ea5"),
                            (0.16, 0.24, "#9bb7d4")):
            ax[0].axvspan(lo, hi, color=col, alpha=0.10)
        ax[0].text(0.06, ax[0].get_ylim()[1] * 0.86,
                   "size-specific\ncrossover",
                   ha="center", fontsize=6.5, color="#3b6ea5")
        ax[0].set_xlabel("simulated $F_{ST}$"); ax[0].set_ylabel("PC $-$ AF AUROC")
        ax[0].set_title("(a) Crossover in one parameter"); ax[0].legend(fontsize=8)

    # (b) governed vs ungoverned
    r = L("s23_governed.json")
    if r:
        a = r["agg"]
        grp = ["AF", "PC-dist"]
        ung = [a["ungoverned_af"][0], a["ungoverned_pc_dist"][0]]
        gov = [a["governed_af"][0], a["governed_pc_dist"][0]]
        x = np.arange(len(grp)); bw = 0.35
        ax[1].bar(x - bw / 2, ung, bw, label="ungoverned", color="#a5523b")
        ax[1].bar(x + bw / 2, gov, bw, label="governed", color="#3b6ea5")
        ax[1].axhline(0.5, ls="--", c="grey", lw=1)
        ax[1].set_xticks(x); ax[1].set_xticklabels(grp)
        ax[1].set_ylabel("AUROC"); ax[1].set_title("(b) Governed egress")
        ax[1].legend(fontsize=8)

    # (c) cross-run composition
    r = L("s24_crossrun.json")
    if r:
        order = r["order"]
        lab = ["$\\lambda_{GC}$", "eig", "AF 1dp", "AF 2dp", "AF full", "loadings"]
        cum = r["agg"]["cumulative"]; mx = r["agg"]["max_single"]
        x = np.arange(len(order))
        ax[2].plot(x, cum, "o-", color="#3b6ea5", label="cumulative")
        ax[2].plot(x, mx, "s--", color="#a5523b", label="max single"
                   )
        ax[2].axhline(0.5, ls=":", c="grey", lw=1)
        ax[2].set_xticks(x); ax[2].set_xticklabels(lab, fontsize=8, rotation=20)
        ax[2].set_ylabel("attack AUROC"); ax[2].set_title("(c) Cross-run composition")
        ax[2].legend(fontsize=8)

    # (d) IBD control
    r = L("s26_ibd.json")
    if r:
        sim = r["sim_control"]["agg"]
        grp = ["IBD planted", "IBD unplanted", "real IBD", "real AF"]
        vals = [sim["planted_ibs"][0], sim["unplanted_ibs"][0],
                r["real"]["chr22"]["agg"]["ibs"][0], r["real"]["chr22"]["agg"]["af"][0]]
        ax[3].bar(range(len(grp)), vals,
                  color=["#3b6ea5", "#9bb7d4", "#a5523b", "#d49b9b"])
        ax[3].axhline(0.5, ls="--", c="grey", lw=1)
        ax[3].set_xticks(range(len(grp)))
        ax[3].set_xticklabels(grp, fontsize=8, rotation=15)
        ax[3].set_ylabel("AUROC"); ax[3].set_ylim(0, 1)
        ax[3].set_title("(d) IBD/haplotype channel")

    fig.savefig(os.path.join(FIG, "fig7_governance.pdf"))
    fig.savefig(os.path.join(FIG, "fig7_governance.png"))
    plt.close(fig)


def fig8_scale():
    """Scale/density/ancestry/IBD panel for Section sec:scale (GPU round)."""
    s27 = merge_shards("s27_gpu_scale_shard*.json")["grid"]
    s35 = load("s35_frontier.json")["grid"]
    s29 = dict(merge_shards("s29_negation_shard*.json"))
    if os.path.exists(os.path.join(RAW, "s29_ladder.json")):
        lad = load("s29_ladder.json")
        for key in ("ladder_sim", "ladder_real", "ladder_msprime"):
            if key in lad:
                s29[key] = lad[key]
    s30 = load("s30_density.json")
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))

    # (a) boundary grid: PC-AF delta vs M, one line per (study size, Fst);
    # Fst is distinguished by line style so no single curve connects different
    # conditions.
    a = ax[0, 0]
    cells = {}
    for k, v in s27.items():
        parts = k.split("_")
        try:
            nsp = int(parts[0].replace("nsp", ""))
            M = int(parts[1].replace("m", ""))
            fst = float(parts[2].replace("f", ""))
        except Exception:
            continue
        cells.setdefault((nsp, fst), []).append((M, v))
    styles = {0.04: "-", 0.08: "--"}
    for (nsp, fst) in sorted(cells):
        rows = sorted(cells[(nsp, fst)], key=lambda r: r[0])
        Ms = [r[0] for r in rows]
        d = [r[1]["delta_pc_minus_af"] for r in rows]
        a.plot(Ms, d, "o" + styles.get(fst, "-"), label=f"$n_{{sp}}$={nsp}, $F_{{ST}}$={fst}")
    a.axhline(0, c="k", lw=0.8)
    a.set_xscale("log")
    a.set_xlabel("markers $M$"); a.set_ylabel("PC $-$ AF AUROC")
    a.set_title("(a) Biobank boundary grid\n(solid: $F_{ST}$=0.04; dashed: $F_{ST}$=0.08)")
    a.legend(fontsize=6, ncol=2)

    # (b) dense small-n frontier
    b = ax[0, 1]
    for nsp in sorted({int(k.split("_")[0].replace("nsp", "")) for k in s35}):
        rows = sorted([(int(k.split("_")[1].replace("m", "")), v)
                       for k, v in s35.items() if k.startswith(f"nsp{nsp}_")])
        b.plot([r[0] for r in rows], [r[1]["delta"] for r in rows], "o-",
               label=f"$n_{{sp}}$={nsp}")
    b.axhline(0, c="k", lw=0.8)
    b.set_xscale("log"); b.set_xlabel("markers $M$"); b.set_ylabel("PC $-$ AF AUROC")
    b.set_title("(b) Dense frontier (small $n$)"); b.legend(fontsize=8)

    # (c) exact k-ladder on sim/real/msprime (log y)
    c = ax[1, 0]
    if "ladder_sim" in s29:
        ks = [k for k in s29["ladder_sim"] if k.isdigit()]
        ks = sorted(ks, key=int)
        c.plot([int(k) for k in ks], [s29["ladder_sim"][k]["mean"] for k in ks],
               "o-", label="sim")
        c.axhline(s29["ladder_sim"]["af_only"]["mean"], ls=":", c="C0", alpha=0.6)
        c.axhline(s29["ladder_sim"]["eig_only"]["mean"], ls="--", c="grey",
                  label="eigen-only (sim)")
    if "ladder_msprime" in s29:
        ks = sorted([k for k in s29["ladder_msprime"] if k.isdigit()], key=int)
        c.plot([int(k) for k in ks], [s29["ladder_msprime"][k]["mean"] for k in ks],
               "s-", label="msprime")
    if "ladder_real" in s29 and isinstance(s29["ladder_real"], dict):
        for chrom, blk in s29["ladder_real"].items():
            if isinstance(blk, dict) and "ladder" in blk:
                ks = sorted(blk["ladder"], key=lambda x: int(x))
                c.plot([int(k) for k in ks], [blk["ladder"][k] for k in ks],
                       "^-", label=chrom)
    c.set_xlabel("top-$k$ PCs"); c.set_ylabel("PC AUROC")
    c.set_title("(c) Exact $k$-ladder"); c.legend(fontsize=7)

    # (d) real density sweep
    d = ax[1, 1]
    for chrom, blk in sorted(s30.items()):
        rows = [(v["n_markers"], v) for v in blk.values() if "n_markers" in v]
        rows.sort()
        d.plot([r[0] for r in rows], [r[1]["pc"] for r in rows], "o-",
               label=f"{chrom} PC")
        d.plot([r[0] for r in rows], [r[1]["af"] for r in rows], "s--",
               label=f"{chrom} AF")
    d.set_xscale("log"); d.set_xlabel("markers (per chr)")
    d.set_ylabel("AUROC"); d.set_title("(d) Real density sweep")
    d.legend(fontsize=7, ncol=2)

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig8_scale.pdf"))
    fig.savefig(os.path.join(FIG, "fig8_scale.png"))
    plt.close(fig)


def merge_shards(pattern):
    """Merge shard files. 'grid'/'paired' are unioned across shards with
    first-wins (setdefault) so a re-run shard cannot flip a cell's value."""
    out = {}
    for f in sorted(glob.glob(os.path.join(RAW, pattern))):
        d = json.load(open(f))
        for key in ("grid", "paired"):
            if key in d:
                tgt = out.setdefault(key, {})
                for k, v in d[key].items():
                    tgt.setdefault(k, v)
        for key in ("ladder_sim", "ladder_real", "ladder_msprime"):
            if key in d:
                out[key] = d[key]
    return out


if __name__ == "__main__":
    fig1_frontier(); fig2_marginal(); fig3_boundary(); fig4_paradox()
    fig5_pc_channel(); fig6_replication(); fig7_governance(); fig8_scale()
    print("figures written to", FIG)
