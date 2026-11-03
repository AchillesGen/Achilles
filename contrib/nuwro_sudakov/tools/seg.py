import sys; sys.path.insert(0,'/home/claude/work/ana')
from common import *
d=load('stepb','')
res={}; zs=[]
for obs,Es,fn in [('sR',[10,30,100,300,900],sigR),('T',[300,1000,2000],trans)]:
    for E in Es:
        row={}
        for g in [0.5,2.0,8.0]:
            t=f'{obs}_b_g{g}_E{E}'
            if t in d:
                v,e=fn(d[t]); a=np.load(f'/home/claude/work/stepb/{t}.npy').astype(float)[:,1]
                row[g]=(v,e,a.mean(),a.std()/np.sqrt(len(a)))
        if len(row)==3:
            vals=np.array([row[g][0] for g in row]); errs=np.array([row[g][1] for g in row]); w=1/errs**2; m=(w*vals).sum()/w.sum(); chi=((vals-m)**2*w).sum()
            n=np.array([row[g][2] for g in row]); ne=np.array([row[g][3] for g in row]); wn=1/ne**2; mn=(wn*n).sum()/wn.sum(); chin=((n-mn)**2*wn).sum()
            res[f'{obs}_{E}']=dict(rows={str(k):v for k,v in row.items()},chi2=chi,chi2_noi=chin)
            print(obs,E,' '.join(f"L={g}: {row[g][0]:.4g}±{row[g][1]:.2g} <N>={row[g][2]:.4f}" for g in row), f"chi2={chi:.2f}/2  chi2(<N>)={chin:.2f}/2")
tot=sum(r['chi2'] for r in res.values()); totn=sum(r['chi2_noi'] for r in res.values())
print('total chi2', tot, 'ndf', 2*len(res), ' <N_int>:', totn)
json.dump(res,open('/home/claude/work/ana/seg.json','w'),indent=1)
