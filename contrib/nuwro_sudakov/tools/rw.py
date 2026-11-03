"""Step 5: reweighting closure.  Nominal sample with per-event weights vs regenerated samples."""
import sys; sys.path.insert(0,'/home/claude/work/ana')
from common import *
VARS=["mfp=0.6","mfp=0.7","mfp=0.82","mfp=1.0","xs=0.8","xs=0.9","xs=1.1","xs=1.25"]
REGEN={"mfp=0.6":"mfp0.6","mfp=1.0":"mfp1.0","xs=0.8":"xs0.8","xs=1.25":"xs1.25"}
def obs(a):
    # a: ana columns: 0 Ek,1 noi,2 nod0,3 nod2,4 nod3,5 pion-int,... 12 npost
    noi=a[:,1]
    return {'P0 (no interaction)':(noi==0)*1.,'<N_int>':noi,'P(N_int>=2)':(noi>=2)*1.,'<N_out>':a[:,12],'<N_elastic>':a[:,2]}
def wmean(x,w):
    m=(w*x).sum()/w.sum()
    # variance of ratio estimator
    n=len(x); v=((w*(x-m))**2).sum()/w.sum()**2
    return m,np.sqrt(v)
res={}
for alg in ['b','a']:
    for cfg in ['T300','T1000','sR50','sR300']:
        nom=f'/home/claude/work/rw/rw_{alg}_{cfg}_nom'
        if not os.path.exists(nom+'.npy') or not os.path.exists(nom+'.w.txt'): continue
        a=np.load(nom+'.npy').astype(float)
        W=np.loadtxt(nom+'.w.txt',comments='#')
        assert len(W)==len(a), (len(W),len(a))
        valid=W[:,1]; W=W[:,2:]
        O=obs(a)
        for iv,v in enumerate(VARS):
            w=W[:,iv]; ess=w.sum()**2/(w**2).sum()
            row=dict(alg=alg,cfg=cfg,var=v,ess_frac=ess/len(w),wmin=float(w.min()),wmax=float(w.max()),wmean=float(w.mean()))
            if v in REGEN:
                g=f'/home/claude/work/rw/rw_{alg}_{cfg}_{REGEN[v]}.npy'
                if os.path.exists(g):
                    b=np.load(g).astype(float); Ob=obs(b)
                    for k in O:
                        m,e=wmean(O[k],w); mb=Ob[k].mean(); eb=Ob[k].std()/np.sqrt(len(b))
                        mn=O[k].mean()
                        row[k]=dict(nominal=mn,reweighted=m,rw_err=e,regenerated=mb,regen_err=eb,z=(m-mb)/np.hypot(e,eb))
            res[f'{alg}_{cfg}_{v}']=row
for k,r in res.items():
    s=f"{k:22s} ESS/N={r['ess_frac']:.3f} w=[{r['wmin']:.3g},{r['wmax']:.3g}] <w>={r['wmean']:.4f}"
    if 'P0 (no interaction)' in r:
        s+=" | "+"; ".join(f"{o}: nom {r[o]['nominal']:.4f} rw {r[o]['reweighted']:.4f}±{r[o]['rw_err']:.4f} regen {r[o]['regenerated']:.4f}±{r[o]['regen_err']:.4f} z={r[o]['z']:+.1f}" for o in ['P0 (no interaction)','<N_int>','<N_out>'])
    print(s)
json.dump(res,open('/home/claude/work/ana/rw.json','w'),indent=1)
