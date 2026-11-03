import json, glob, os, numpy as np
R_FB = 8.00043364            # fm, fb_6_6a_pointDensity[0]  [read in source, nucleus_data.cc:232]
AREA_MB = np.pi*(0.999*R_FB)**2*10.0   # beam disk of radius 0.999 R (nucleusmaker.h:41), 1 fm^2 = 10 mb
def load(dirn, prefix):
    out={}
    for f in glob.glob(f'/home/claude/work/{dirn}/{prefix}*.sum.json'):
        s=json.load(open(f)); out[s['tag']]=s
    return out
def pint(s):
    n=s['nev']; k=s['nint']; p=k/n; return p, np.sqrt(p*(1-p)/n)
def sigR(s):
    p,e=pint(s); return AREA_MB*p, AREA_MB*e
def trans(s):
    p,e=pint(s); return 1-p, e
def wfit(x,y,e,deg):
    """weighted polynomial fit; returns coeffs (c0..), covariance, chi2, ndf"""
    x=np.asarray(x); y=np.asarray(y); e=np.asarray(e)
    A=np.vander(x,deg+1,increasing=True)/e[:,None]; b=y/e
    cov=np.linalg.inv(A.T@A); c=cov@A.T@b
    chi2=float(((A@c-b)**2).sum()); return c,cov,chi2,len(x)-deg-1
