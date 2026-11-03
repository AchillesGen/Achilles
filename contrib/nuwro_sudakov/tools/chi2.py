import sys; sys.path.insert(0,'/home/claude/work/ana')
from common import *
import csv
def readcsv(f):
    rows=[r for r in csv.reader(l for l in open(f) if not l.startswith('#'))]
    return np.array([[float(x) for x in r] for r in rows[1:]])
dR=readcsv('/home/claude/work/data/pC_reaction_xsec_Carlson.csv')
dT=readcsv('/home/claude/work/data/C12_transparency.csv')
def model_curve(dirn, prefix, fn):
    d=load(dirn,prefix); pts=sorted((s['Ek'],)+fn(s) for s in d.values())
    return np.array(pts)
out={}
for v in ['a','b']:
    mR=model_curve('chi',f'chi_sR_{v}_',sigR); mT=model_curve('chi',f'chi_T_{v}_',trans)
    if len(mR)<5: continue
    # sigma_R: interpolate in log E (model grid covers 10..860; data down to 9.9 -> clamp)
    E=np.clip(dR[:,0],mR[0,0],mR[-1,0])
    yR=np.interp(np.log(E),np.log(mR[:,0]),mR[:,1]); eR=np.interp(np.log(E),np.log(mR[:,0]),mR[:,2])
    pulls=(yR-dR[:,1])/dR[:,2]
    res=dict(sR_chi2=float((pulls**2).sum()),sR_N=len(pulls))
    for lo,hi in [(0,100),(100,1000)]:
        m=(dR[:,0]>=lo)&(dR[:,0]<hi); res[f'sR_chi2_{lo}_{hi}']=float((pulls[m]**2).sum()); res[f'sR_N_{lo}_{hi}']=int(m.sum())
        res[f'sR_ratio_{lo}_{hi}']=float(np.mean(yR[m]/dR[m,1]))
    if len(mT)>=5:
        yT=np.interp(dT[:,0],mT[:,0],mT[:,1]); pT=(yT-dT[:,1])/dT[:,2]
        res.update(T_chi2=float((pT**2).sum()),T_N=len(pT),T_ratio=float(np.mean(yT/dT[:,1])))
    res['model_stat_rel_sR']=float(np.mean(eR/yR))
    out[v]=res; print(v,res)
json.dump(out,open('/home/claude/work/ana/chi2.json','w'),indent=1)
