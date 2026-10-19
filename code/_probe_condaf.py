import sys, numpy as np
sys.path.insert(0,'/data/lvvyongji1/GWAS/code')
from common import (simulate_population, split_population3, release_objects,
                    allele_freq, af_lrt_attack, auroc)
for seed in (0,1,2):
    pop=simulate_population(n_groups=3,n_per_group=1200,n_snps=4000,seed=seed)
    site,ref,non=split_population3(pop,300,300,300)
    rel=release_objects(site); mask=rel["mask"]; af_s=rel["af"]; af_r=allele_freq(ref["G"])
    G=site["G"]; N=G.shape[0]; anc=site["ancestry"]
    rng=np.random.default_rng(seed+999); ts=rng.random(N)<0.5
    K=np.where(ts)[0]; T=np.where(~ts)[0]
    def a_of(afvec, idx):
        am=af_lrt_attack(G,afvec,af_r,mask=mask); an=af_lrt_attack(non["G"],afvec,af_r,mask=mask)
        y=np.concatenate([np.ones(len(idx)),np.zeros(len(non["G"]))])
        return float(auroc(np.concatenate([am[idx],an]),y))
    def residual(remove):
        num=G.sum(0).astype(float)-G[remove].sum(0).astype(float)
        return np.clip(num/(2.0*(N-len(remove))),1e-6,1-1e-6)
    # V1 random-removal residual, eval test
    print(seed,"V1 rand-remove resid, eval test  :", round(a_of(residual(K),T),4))
    # V2 balanced-per-ancestry removal (half of each ancestry), eval ALL members
    half=[];
    for g in np.unique(anc):
        gi=np.where(anc==g)[0]; half.append(gi[:len(gi)//2])
    half=np.concatenate(half)
    print(seed,"V2 balanced half-remove, eval all :", round(a_of(residual(half),np.arange(N)),4))
    print(seed,"V2b balanced half-remove, eval rest:", round(a_of(residual(half),np.setdiff1d(np.arange(N),half)),4))
    # V3 plain AF eval test (baseline)
    print(seed,"V3 plain AF, eval test            :", round(a_of(af_s,T),4))
    # V4 plain AF eval all
    print(seed,"V4 plain AF, eval all             :", round(a_of(af_s,np.arange(N)),4))
