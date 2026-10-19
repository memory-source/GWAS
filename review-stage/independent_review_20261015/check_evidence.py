"""Bounded diagnostics for independent review. No full study experiments are rerun."""
import ast, hashlib, json, os, re, subprocess, sys, tempfile
from pathlib import Path
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

def load(name):
    return json.loads((ROOT / "results/raw" / (name + ".json")).read_text())

def hook_probe():
    source = (ROOT / "code/enforce.py").read_text()
    names = {"SandboxViolation", "install_egress_audit_hook"}
    hook_source = "\n\n".join(ast.get_source_segment(source, n) for n in ast.parse(source).body if getattr(n, "name", None) in names)
    program = """
import os, sys, json
import numpy as np
WRITE_BITS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND
exec(os.environ["REVIEW_HOOK_SOURCE"])
fd = int(os.environ["REVIEW_REPORT_FD"])
base = os.environ["REVIEW_TEMP_DIR"]
install_egress_audit_hook([fd])
report = {}
try:
    with open(os.path.join(base, "ordinary_open.bin"), "wb") as f:
        f.write(b"SYNTHETIC_ONLY")
    report["ordinary_open"] = "allowed"
except Exception as e:
    report["ordinary_open"] = type(e).__name__
try:
    np.array([41.0, 42.0], dtype=np.float64).tofile(os.path.join(base, "numpy_tofile.bin"))
    report["numpy_tofile"] = "allowed"
except Exception as e:
    report["numpy_tofile"] = type(e).__name__
try:
    os.mkdir(os.path.join(base, "SYNTHETIC_VALUE_042"))
    report["ordinary_mkdir"] = "allowed"
except Exception as e:
    report["ordinary_mkdir"] = type(e).__name__
os.write(fd, json.dumps(report).encode())
os.close(fd)
"""
    with tempfile.TemporaryDirectory(prefix="gwas_review_synthetic_") as temp:
        read_fd, write_fd = os.pipe()
        env = dict(os.environ, REVIEW_HOOK_SOURCE=hook_source, REVIEW_REPORT_FD=str(write_fd), REVIEW_TEMP_DIR=temp, PYTHONDONTWRITEBYTECODE="1")
        proc = subprocess.Popen([sys.executable, "-c", program], env=env, pass_fds=(write_fd,), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.close(write_fd)
        with os.fdopen(read_fd, "rb") as stream:
            raw = stream.read()
        rc = proc.wait(timeout=30)
        result = json.loads(raw) if raw else {}
        p = Path(temp) / "numpy_tofile.bin"
        result.update(child_exit=rc, ordinary_open_file_exists=(Path(temp) / "ordinary_open.bin").exists(), numpy_tofile_file_exists=p.exists(), numpy_tofile_bytes=p.stat().st_size if p.exists() else 0, synthetic_directory_exists=(Path(temp) / "SYNTHETIC_VALUE_042").is_dir(), scope="Exact hook in a fresh child; not a full study/runtime replay. Synthetic temporary data only.")
        return result

def main():
    tex = (ROOT / "paper/main.tex").read_text()
    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S).group(1)
    pairs = load("s10_paired")
    rows = pairs["per_seed"]
    means = {k: float(np.mean([r[k] for r in rows])) for k in ("af_only", "pc_calibrated", "af_conditioned", "combined")}
    cohorts = load("s48_real_structure")["cohorts"]
    small = {n: {k: v[k] for k in ("n_study", "M", "af_auc", "pc_auc", "gap")} for n,v in cohorts.items() if n.startswith("pop:")}
    largest = {str(chrom): {k:v for k,v in load("s56_real_poolmax_chr"+str(chrom))["cells"].items() if "350" in k} for chrom in (22,1,2)}
    z = np.load(ROOT / "results/cache/real_panel_chr22_c284b941.npz")
    meta = [line.split() for line in (ROOT / "data/external/1000g_panel.txt").read_text().splitlines()[1:] if line.strip()]
    pops = np.array([r[1] for r in meta])
    supers = np.array([r[2] for r in meta])
    if not np.array_equal(z["sup"], supers):
        raise ValueError("Cache and panel labels are not aligned")
    definition_check = dict(cache_shape=list(z["G"].shape), observed_position_min=int(z["pos"].min()), observed_position_max=int(z["pos"].max()), scope="Checks pooled-label identity and original label counts; not an estimate that genetic structure is absent.", superpopulations={})
    for sp in ("EUR", "EAS", "AFR"):
        take = supers == sp
        labels, counts = np.unique(pops[take], return_counts=True)
        af = z["G"][take].mean(0) / 2
        ht = 2 * af * (1 - af)
        definition_check["superpopulations"][sp] = dict(original_population_counts={str(k):int(v) for k,v in zip(labels,counts)}, fst_like_after_single_label_pooling=float(np.mean((ht[ht > 0]-ht[ht > 0])/ht[ht > 0])))
    d = np.array([r["combined"]-r["af_only"] for r in rows])
    rng = np.random.default_rng(0)
    boot = [float(d[rng.integers(0,len(d),len(d))].mean()) for _ in range(5000)]
    recomputed_ci = [float(x) for x in np.percentile(boot,[2.5,97.5])]
    sys.path.insert(0, str(ROOT / "code"))
    from common import auroc
    paths = ("paper/main.tex","paper/main.pdf","code/enforce.py","code/common.py","code/s25_paradox.py","results/raw/s10_paired.json")
    output = dict(scope="Archived-value recomputation and bounded hook probes; no full experiment reproduction.", inputs_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}, abstract_whitespace_word_count=len(abstract.split()), paired_s10=dict(n_seeds=len(rows), recomputed_means=means, max_abs_mean_error_against_archive=max(abs(means[k]-pairs["means"][k]) for k in means), combined_minus_unconditioned_af_recomputed=means["combined"]-means["af_only"], combined_minus_unconditioned_af_archived_ci=pairs["combined_vs_af"]["gap_ci95"], conditioned_af_beats_pc_seed_count=sum(r["af_conditioned"]>r["pc_calibrated"] for r in rows), combined_minus_conditioned_af_mean=means["combined"]-means["af_conditioned"]), self_reference_bias={k:v for k,v in load("s1d_bias").items() if k!="rows"}, real_single_population_cohorts=small, real_largest_cells=largest, cpu_constant_score_auc=float(auroc(np.ones(4),np.array([1,1,0,0]))), governance_archived_agg=load("s39_enforcement")["agg"], audit_hook_probe=hook_probe())
    import pymupdf
    with pymupdf.open(ROOT/"paper/main.pdf") as pdf:
        output.update(pdf_pages=len(pdf), pdf_creation_date=pdf.metadata.get("creationDate"))
    output["structure_definition_check"] = definition_check
    output["paired_s10"]["combined_minus_unconditioned_af_recomputed_ci"] = recomputed_ci
    output["paired_s10"]["ci_max_abs_error_against_archive"] = max(abs(a-b) for a,b in zip(recomputed_ci,pairs["combined_vs_af"]["gap_ci95"]))
    path = OUT/"evidence_checks.json"
    path.write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({k:output[k] for k in ("paired_s10","audit_hook_probe","pdf_pages")},ensure_ascii=False,indent=2))
if __name__ == "__main__":
    main()
