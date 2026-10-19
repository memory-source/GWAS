"""Read-only replay of the EAS random-span control; writes only this review evidence."""
from pathlib import Path
import hashlib,json,os,sys
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
sys.dont_write_bytecode=True
import numpy as np
ROOT=Path("/data/lvvyongji1/GWAS")
sys.path.insert(0,str(ROOT/"code"))
import s61_spectral_tomography as s
cache=ROOT/"results/cache/real_panel_chr22_e6f1c0cf.npz"
with np.load(cache,allow_pickle=True) as z:
    G=z["G"][z["sup"]=="EAS"]
groups=np.repeat("EAS",len(G))
rows=[]
role_sets=[]
for seed in range(5):
    roles=s.split_roles(groups,s.CAPS,seed)
    honest,partner,nonmem=[G[roles[t]] for t in ("honest","partner","nonmember")]
    mu,sd=s.joint_baseline(honest,partner)
    y=np.r_[np.ones(len(honest)),np.zeros(len(nonmem))]
    cal=np.random.default_rng(seed+31).random(len(y))<0.5
    R=np.random.default_rng(seed+71).standard_normal((G.shape[1],10))
    unclipped=s._membership_auc_from_operator(R,honest,nonmem,mu,sd,y,cal,z_clip=None)[0]
    clipped=s._membership_auc_from_operator(R,honest,nonmem,mu,sd,y,cal,z_clip=5.0)[0]
    rows.append(dict(seed=seed,unclipped=unclipped,clipped=clipped,test_nonmembers=int((~cal[len(honest):]).sum())))
    role_sets.append(set(roles["honest"]))
out=dict(scope="Same EAS cache, splits and random spans; only z_clip changes. This tests control comparability, not invalidity of the active-versus-passive gap.",cache=str(cache.relative_to(ROOT)),cache_sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),EAS_shape=list(G.shape),per_seed=rows,mean_unclipped=float(np.mean([r["unclipped"] for r in rows])),mean_clipped=float(np.mean([r["clipped"] for r in rows])),honest_overlap_seed0_others=[len(role_sets[0]&t) for t in role_sets[1:]])
outpath=Path(__file__).resolve().parent/"control_clipping_checks.json"
outpath.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
print(json.dumps(out,ensure_ascii=False,indent=2))
