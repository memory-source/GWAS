"""Bounded independent archive checks and synthetic current-hook probe.
Research artifacts are read-only; this script writes evidence_checks.json beside itself.
"""
import ast, hashlib, json, os, re, subprocess, sys, tempfile
from pathlib import Path
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
def load(name):
    return json.loads((ROOT / "results/raw" / (name + ".json")).read_text())
def probe():
    source = (ROOT / "code/enforce.py").read_text()
    extracted = "\n\n".join(ast.get_source_segment(source,n) for n in ast.parse(source).body if getattr(n,"name",None) in {"SandboxViolation","install_egress_audit_hook"})
    program = """
import os, sys, json
WRITE_BITS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND
exec(os.environ["REVIEW_HOOK"])
fd = int(os.environ["REVIEW_FD"])
base = os.environ["REVIEW_TEMP"]
install_egress_audit_hook([fd])
report = {}
try:
    with open(os.path.join(base,"ordinary_write.bin"),"wb") as f:
        f.write(b"SYNTHETIC_ONLY")
    report["ordinary_write"] = "allowed"
except Exception as exc:
    report["ordinary_write"] = type(exc).__name__
try:
    os.mkdir(os.path.join(base,"SYNTHETIC_VALUE_042"))
    report["ordinary_mkdir"] = "allowed"
except Exception as exc:
    report["ordinary_mkdir"] = type(exc).__name__
os.write(fd,json.dumps(report).encode())
os.close(fd)
"""
    with tempfile.TemporaryDirectory(prefix="gwas_review_20261017_") as tmp:
        rfd,wfd=os.pipe()
        env=dict(os.environ,REVIEW_HOOK=extracted,REVIEW_FD=str(wfd),REVIEW_TEMP=tmp,PYTHONDONTWRITEBYTECODE="1")
        p=subprocess.Popen([sys.executable,"-c",program],env=env,pass_fds=(wfd,),stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        os.close(wfd)
        try:
            _,err=p.communicate(timeout=20)
        except subprocess.TimeoutExpired:
            p.kill()
            p.communicate()
            raise
        with os.fdopen(rfd,"rb") as f:
            raw=f.read()
        if p.returncode or not raw:
            raise RuntimeError(f"probe failed: {p.returncode}; {err.decode()}")
        result=json.loads(raw)
        result.update(child_exit=p.returncode,ordinary_file_exists=(Path(tmp)/"ordinary_write.bin").exists(),encoded_directory_exists=(Path(tmp)/"SYNTHETIC_VALUE_042").is_dir(),scope="Exact current hook in fresh subprocess; synthetic temporary data; no native-code loading; not full runtime replay.")
        return result
def synthetic_algebra():
    import numpy as np
    rng=np.random.default_rng(117)
    honest=rng.integers(0,3,size=(40,24)).astype(float)
    own=rng.integers(0,3,size=(20,24)).astype(float)
    mu=np.vstack([honest,own]).mean(0)
    H=honest-mu
    CH=H.T@H
    max_residual=0.0
    hist_preserved=True
    Qs,Ys=[],[]
    ranks=[]
    for t in range(8):
        partner=own.copy()
        if t:
            a,b=rng.choice(len(own),2,replace=False)
            block=int(rng.integers(0,6))
            sl=slice(4*block,4*(block+1))
            partner[a,sl],partner[b,sl]=own[b,sl],own[a,sl]
        hist_preserved=hist_preserved and all(np.array_equal((own==j).sum(0),(partner==j).sum(0)) for j in range(3))
        P=partner-mu
        K=P.T@P
        lam,V=np.linalg.eigh(CH+K)
        V=V[:,-3:]
        lam=lam[-3:]
        Y=V*lam[None,:]-K@V
        max_residual=max(max_residual,float(np.linalg.norm(Y-CH@V)/np.linalg.norm(CH@V)))
        Qs.append(V)
        Ys.append(Y)
        ranks.append(int(np.linalg.matrix_rank(np.concatenate(Qs,axis=1))))
    Q=np.concatenate(Qs,axis=1)
    Y=np.concatenate(Ys,axis=1)
    return dict(synthetic_only=True,identity_max_relative_residual=max_residual,histograms_preserved=hist_preserved,stacked_ranks=ranks,final_scatter_relative_error=float(np.linalg.norm(Y@np.linalg.pinv(Q)-CH)/np.linalg.norm(CH)))

def main():
    import numpy as np
    pair=load("s10_paired")
    rows=pair["per_seed"]
    chans=("af_only","pc_calibrated","af_conditioned","combined")
    means={c:float(np.mean([r[c] for r in rows])) for c in chans}
    d=np.array([r["combined"]-r["af_only"] for r in rows])
    rng=np.random.default_rng(0)
    ci=np.percentile([d[rng.integers(0,len(d),len(d))].mean() for _ in range(5000)],[2.5,97.5]).tolist()
    tex=(ROOT/"paper/main.tex").read_text()
    abstract=re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}",tex,re.S).group(1)
    plink=load("s63_plink_export_chr22")
    vals=[v for cell in load("s66_control_provenance_chr22")["cohorts"].values() for k,v in cell.items() if not k.endswith("_ci95")]
    paths=["paper/main.tex","code/enforce.py","code/s10_paired.py","code/s61_spectral_tomography.py","code/s63_plink_export.py","code/s65_identity_ablation.py","code/s66_control_provenance.py","results/raw/s10_paired.json","results/raw/s61_spectral_tomography_chr22.json","results/raw/s61_spectral_tomography_chr1.json","results/raw/s63_plink_export_chr22.json"]
    out=dict(scope="Archive recomputation, code inspection and synthetic hook probe; not full independent PCA/GWAS/PLINK reproduction.",
        inputs_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        paired_s10=dict(n_seeds=len(rows),recomputed_means=means,max_mean_error_against_archive=max(abs(means[c]-pair["means"][c]) for c in chans),conditioned_af_beats_pc_seeds=sum(r["af_conditioned"]>r["pc_calibrated"] for r in rows),combined_minus_conditioned_af=means["combined"]-means["af_conditioned"],combined_minus_af_recomputed_ci95=ci,ci_max_error_against_archive=max(abs(a-b) for a,b in zip(ci,pair["combined_vs_af"]["gap_ci95"]))),
        plink_export_archived={k:plink[k] for k in ("plink_version","normalization","identity","gap_q8_mean","gap_q8_ci95","gap_q8_pos_splits")},
        raw_v_ablation_archived=load("s65_identity_ablation_chr22")["cohorts"],
        no_release_controls_archived=dict(n_means=len(vals),min_mean=min(vals),max_mean=max(vals)),
        manuscript=dict(abstract_whitespace_words=len(abstract.split()),contains_active_branch_names=any(s in tex for s in ("s61_","s62_","s63_","tomograph"))),
        synthetic_hook_probe=probe(),synthetic_algebra=synthetic_algebra())
    (OUT/"evidence_checks.json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({k:out[k] for k in ("paired_s10","no_release_controls_archived","manuscript","synthetic_hook_probe","synthetic_algebra")},indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
