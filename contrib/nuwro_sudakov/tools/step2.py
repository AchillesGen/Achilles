import sys; sys.path.insert(0,'/home/claude/work/ana')
from common import *
d=load('step2','')
steps=[0.00625,0.0125,0.025,0.05,0.1,0.2,0.4,0.8,1.6]
res={}
for obs,Es,fn in [('sR',[10,20,30,50,70,100,150,200,300,400,550,700,900],sigR),('T',[200,300,450,600,800,1000,1300,1600,2000],trans)]:
    for E in Es:
        xs=[];ys=[];es=[]
        for st in steps:
            t=f"{obs}_a_s{st}_E{E}"
            if t in d: v,e=fn(d[t]); xs.append(st); ys.append(v); es.append(e)
        if len(xs)<4: continue
        xs=np.array(xs); ys=np.array(ys); es=np.array(es)
        # quadratic fit over all steps, linear fit over steps <= 0.4
        mq=xs<=0.8
        c2,cov2,chi2q,ndfq=wfit(xs[mq],ys[mq],es[mq],2)
        m=xs<=0.4
        c1,cov1,chi2l,ndfl=wfit(xs[m],ys[m],es[m],1) if m.sum()>=3 else (np.array([np.nan,np.nan]),np.full((2,2),np.nan),np.nan,0)
        g=np.array([0,0.2,0.04]); bias=c2@g; bias_e=np.sqrt(g@cov2@g)
        res[f"{obs}_{E}"]=dict(obs=obs,E=E,steps=xs.tolist(),vals=ys.tolist(),errs=es.tolist(),
            zero_q=c2[0],zero_q_e=np.sqrt(cov2[0,0]),chi2q=chi2q,ndfq=ndfq,
            zero_l=c1[0],zero_l_e=np.sqrt(cov1[0,0]),chi2l=chi2l,ndfl=ndfl,
            slope_l=c1[1],slope_l_e=np.sqrt(cov1[1,1]),
            bias02=bias,bias02_e=bias_e,val02=ys[list(xs).index(0.2)] if 0.2 in xs else None)
        r=res[f"{obs}_{E}"]
        print(f"{obs:3s} E={E:5d}  zero(quad)={r['zero_q']:.4g}±{r['zero_q_e']:.2g} chi2={chi2q:.1f}/{ndfq}  zero(lin<=0.4)={r['zero_l']:.4g}±{r['zero_l_e']:.2g}  bias@0.2={bias:+.3g}±{bias_e:.2g} ({100*bias/r['zero_q']:+.2f}%)")
json.dump(res,open('/home/claude/work/ana/step2.json','w'),indent=1)
