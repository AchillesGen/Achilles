import numpy as np, sys, glob, os, json
from scipy import stats
R=1.2*12**(1/3); rho=1/(4/3*np.pi*1.2**3)
out={}
for d in sorted(glob.glob('/home/claude/work/unit/*/sud_test.txt')):
    tag=d.split('/')[-2]; sig_mb=float(tag.split('_s')[-1]) if '_s1' in tag and tag.endswith(('10','100')) else 40.
    if tag.endswith('_s10'): sig_mb=10.
    if tag.endswith('_s100'): sig_mb=100.
    G=rho*sig_mb/10.  # fm^-1 (1 mb = 0.1 fm^2)
    a=np.loadtxt(d)
    x0,y0,z0,c,xc,yc,zc=a.T
    b2=x0**2+y0**2; L=np.sqrt(R**2-b2)-z0
    p=1-np.exp(-G*L)
    n=len(a); k=c.sum(); ek=p.sum(); sd=np.sqrt((p*(1-p)).sum())
    m=c==1; dpt=zc[m]-z0[m]
    u=(1-np.exp(-G*dpt))/p[m]
    ks=stats.kstest(u,'uniform')
    out[tag]=dict(sigma_mb=sig_mb,n=n,Pint=k/n,Pint_exact=ek/n,z=(k-ek)/sd,KS_D=ks.statistic,KS_p=ks.pvalue,
                  maxdev_xy=float(np.max(np.abs(np.r_[xc[m]-x0[m],yc[m]-y0[m]]))))
    print(f"{tag:22s} sig={sig_mb:5.0f} N={n} P_int={k/n:.5f} exact={ek/n:.5f} z={(k-ek)/sd:+.2f} KS D={ks.statistic:.4f} p={ks.pvalue:.3f}")
json.dump(out,open('/home/claude/work/unit/ks_results.json','w'),indent=1)
