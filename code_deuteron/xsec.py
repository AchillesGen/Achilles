import xsec_fact
import numpy as np
from scipy import interpolate
from mpi4py import MPI

nev = 80000
neq = 10000
nvoid = 10

mp = 938.272046
mn = 939.56563
mu = 931.4940610
hbarc = 197.327053

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    nwlk = 8
    ee = 1511.0
    thetalept = 89.97
    nZ = 1
    nA = 2
    fg = 0
    kf = 225
    nw = 400
    wmin = 800
    wmax = 1200
else:
    nwlk = None
    ee = None
    thetalept = None
    nZ = None
    nA = None
    fg = None
    nw = None
    wmax = None

nwlk = comm.bcast(nwlk, root=0)
ee = comm.bcast(ee, root=0)
thetalept = comm.bcast(thetalept, root=0)
nZ = comm.bcast(nZ, root=0)
nA = comm.bcast(nA, root=0)
fg = comm.bcast(fg, root=0)
nw = comm.bcast(nw, root=0)
wmax = comm.bcast(wmax, root=0)

irn0 = np.array([19+i for i in range(nwlk)])
irn = np.array(irn0[0:nwlk])

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

# Construct omega grid and the form factors
w = np.empty(nw)
wt = np.empty(nw)
Q2 = np.empty(nw)
hw = (wmax-wmin)/nw
for i in range(nw):
    w[i] = i*hw+wmin
    eef = ee - w[i]
    Q2[i] = 2.0*ee*eef*(1.0 - coste)

# Compute the response functions in the impulse approximation
qval = np.sqrt(Q2 + w**2)
for iw in range(nw):
    r_avg = 0.0
    r_err = 0.0
    i_acc = 0.0
    i_avg = 0
    ip_o = np.zeros(nwlk, dtype=int)
    ip_n = np.zeros(nwlk, dtype=int)
    g_o = np.zeros(nwlk)
    g_n = np.zeros(nwlk)
    f_o = np.zeros(nwlk)

    for i in range(nwlk):
        while g_o[i] <= 0.0:
            ip_o[i] = int(n_p*np.random.random()-1)
            g_o[i] = g_eval(p[ip_o[i]],pke[0,ip_o[i]])

    for iv in range(nev):
        for j in range(nwlk):
            ip_n[j] = int(ip_o[j]+0.5*n_p*(-1+2*np.random.random())-1)
            if ip_n[j] < n_p and ip_n[j] >= 0 :
                g_n[j] = g_eval(p[ip_n[j]], pke[0,ip_n[j]])
            else:
                g_n[j] = 0.0
            if g_n[j]/g_o[j] >= np.random.rand():
                ip_o[j] = ip_n[j]
                g_o[j] = g_n[j]
                i_acc += 1
#            print(ip_o[j],ip_n[j])

            if iv >= neq and iv%nvoid == 0:
                f_o[j] = f_eval(nZ,p[ip_o[j]],w[iw],qval[iw],mqe,pke[0,ip_o[j]],fg,kf,
                                thetalept,ee,iform)
                f_o[j] *= 1e9/g_o[j]
                r_avg += f_o[j]
                r_err += f_o[j]**2
                i_avg += 1

    r_avg /= i_avg
    r_err /= i_avg
    r_err = np.sqrt(r_err-r_avg**2)/(i_avg-1)
    print(w[iw], r_avg, r_err)





