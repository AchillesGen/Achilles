import xsec_fact
import numpy as np

mp = 938.272046
mn = 939.56563
hbarc = 197.327053

mqe = 0.5*(mp+mn)
mqef = mqe/hbarc
xsec_fact.dirac_matrices.dirac_matrices_in(mqef)

def f_eval(nZ,p,w,qval,mqe,pke,fg,kf,thetalept,ee,iform):
    ep = np.sqrt(p**2+mqe**2)
    if fg == 1:
        wt = w
    else:
        wt = w-2

    cost_te = ((wt+ep)**2 - p**2 - qval**2 - mqe**2)/(2*p*qval)
    if abs(cost_te) > 1:
        return 0

    pf = np.sqrt(p**2 + qval**2 + 2*qval*p*cost_te)
    epf = np.sqrt(mqe**2 + pf**2)
    if pf > kf:
        phi = 2.0*np.pi*np.random.rand(1)
        sig = xsec_fact.cc1(qval/hbarc,w,wt,p/hbarc,pf/hbarc,phi,ee,thetalept,iform) 
        return p**2*pke*(nZ*sig)*epf/(p*qval)*2*np.pi
    return 0

outpy = f_eval(1,472.85317271179161,100,908.0056717,938.918838,
               2.3725323720921285e-11,0,225,0.13962634015954636,6519,1)

outf90 = 6.8167880495597444E-011 

print('Fortran: {}\nPython: {}\n|1-Python/Fortran| = {}'.format(outf90,outpy,abs(1-outf90/outpy)))
