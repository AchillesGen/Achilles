import sys; sys.path.insert(0,'/home/claude/work/ana')
from common import *
steps=[0.00625,0.0125,0.025,0.05,0.1,0.2,0.4,0.8,1.6]
res={}
for obs,Es in [('sR',[10,30,50,100,200,300,550,900]),('T',[200,300,600,1000,2000])]:
    for E in Es:
        xs=[];ys=[];es=[]
        for st in steps:
            f=f'/home/claude/work/step2/{obs}_a_s{st}_E{E}.npy'
            if os.path.exists(f):
                a=np.load(f).astype(float); x=a[:,1]; xs.append(st); ys.append(x.mean()); es.append(x.std()/np.sqrt(len(x)))
        xs=np.array(xs); ys=np.array(ys); es=np.array(es)
        m=xs<=0.4; c,cov,chi2,ndf=wfit(xs[m],ys[m],es[m],2 if m.sum()>=4 else 1)
        fb=f'/home/claude/work/stepb/{obs}_b_g2.0_E{E}.npy'
        b=np.load(fb).astype(float)[:,1]; vb=b.mean(); eb=b.std()/np.sqrt(len(b))
        v02=ys[list(xs).index(0.2)]
        res[f'{obs}_{E}']=dict(steps=xs.tolist(),vals=ys.tolist(),errs=es.tolist(),zero=c[0],zero_e=np.sqrt(cov[0,0]),b=vb,b_e=eb,a02=v02)
        print(f"{obs} {E:5d}: <N_int> a(0.2)={v02:.4f}  a(0)={c[0]:.4f}±{np.sqrt(cov[0,0]):.4f} (chi2 {chi2:.1f}/{ndf})  b={vb:.4f}±{eb:.4f}  a0.2 bias={100*(v02/c[0]-1):+.1f}%  b-a0 z={(vb-c[0])/np.hypot(eb,np.sqrt(cov[0,0])):+.1f}")
json.dump(res,open('/home/claude/work/ana/noi.json','w'),indent=1)
