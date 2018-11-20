import xsec_fact
import numpy as np
from scipy import interpolate

mp = 938.272046
mn = 939.56563
mu = 931.4940610
hbarc = 197.327053

nwlk = 1
ee = 1511.0
thetalept = 89.97
nZ = 1
nA = 2
fg = 0
nw = 400
wmax = 1200

irn0 = [19+i for i in range(nwlk)]
irn = irn0[0:nwlk]

iform = 1
mqe = 0.5*(mp+mn)
mqef = mqe/hbarc
mnuc = nA*mu
pi = np.pi
thetalept=thetalept*pi/180.0
coste = np.cos(thetalept)
xsec_fact.dirac_matrices.dirac_matrices_in(mqef)

def g_eval(pj,pke):
    return (4.0*np.pi)*pj**2*pke 

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

if fg != 1:
    with open('nv2IIB.dat') as f:
        nskip = 2
        nbox = 201
        p0 = np.empty(nbox)
        pke0 = np.empty((1,nbox))
        for i, line in enumerate(f):
            if i < nskip:
                continue
            tokens = line.split()
            p0[i-nskip] = float(tokens[0])*hbarc
            pke0[0,i-nskip] = float(tokens[-1])/hbarc**3/(2*pi)**3

    pmax = p0[-1]
    n_p = nbox*6
    ne = 1
    p = np.empty(n_p)
    pke = np.empty((ne,n_p))
    dp = np.empty(n_p)
    xe = np.empty(ne)
    xe[0] = 2.0
    hp = pmax/n_p
    for i in range(n_p):
        p[i] = (i+0.5)*hp

    p_interp = interpolate.interp1d(p0,pke0[0],kind='cubic')
    pke[0] = p_interp(p)

    for i in range(n_p):
        dp[i] = sum(pke[:,i])
    norm = sum(p**2*dp)*4*pi*hp
    pke /= norm
    print('n(k) norm initial = {0}'.format(norm))

    for i in range(n_p):
        dp[i] = sum(pke[:,i])
    norm = sum(p**2*dp)*4*pi*hp
    print('n(k) norm = {0}'.format(norm))

    print(f_eval(1,472.85317271179161,100,908.0056717,938.918838,
        2.3725323720921285e-11,0,225,0.13962634015954636,6519,1))

# Construct omega grid and the form factors
    w = np.empty(nw)
    wt = np.empty(nw)
    Q2 = np.empty(nw)
    hw = wmax/nw
    for i in range(nw):
        w[i] = (i+1)*hw
        eef = ee - w[i]
        Q2[i] = 2.0*ee*eef*(1.0 - coste)

# Compute the response functions in the impulse approximation
    qval = np.sqrt(Q2 + w**2)
    for iw in range(nw):
        r_avg = 0.0
        r_err = 0.0
        i_acc = 0.0
        i_avg = 0
        g_o = 0.0

        for i in range(nwlk):
            while g_o <= 0.0:
                ip_o = int(n_p*np.random.rand(1))
                g_o = g_eval(p[ip_o],pke[0,ip_o])



