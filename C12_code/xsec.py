import xsec_fact
import numpy as np
from scipy import interpolate
import vegas
import matplotlib.pyplot as plt

nev = 80000
neq = 10000
nvoid = 10

mp = 938.272046
mn = 939.56563
mu = 931.4940610
hbarc = 197.327053

nwlk = 8
ee = 730.0
thetalept = 37.1
nZ = 6
nA = 12
fg = 0
kf = 225
nw = 20
wmin = 5
wmax = 650

iform = 2
mqe = 0.5*(mp+mn)
mqef = mqe/hbarc
thetalept=thetalept*np.pi/180.0
coste = np.cos(thetalept)
xsec_fact.dirac_matrices.dirac_matrices_in(mqef)

def g_eval(pj,pke):
    return (4.0*np.pi)*pj**2*pke 

def f_eval(nZ,p,e,w,qval,mqe,pke,fg, kf,thetalept,ee,iform):
    ep = np.sqrt(p**2+mqe**2)
    if fg == 1:
        wt = w
    else:
        wt = w-e+mqe-ep

    cost_te = ((wt+ep)**2 - p**2 - qval**2 - mqe**2)/(2*p*qval)
    if abs(cost_te) > 1:
        return 0

    pf = np.sqrt(p**2 + qval**2 + 2*qval*p*cost_te)
    epf = np.sqrt(mqe**2 + pf**2)
    if pf >= kf:
        phi = 2.0*np.pi*np.random.rand(1)
        sig = xsec_fact.cc1(qval/hbarc,w,wt,p/hbarc,pf/hbarc,phi,ee,thetalept,iform) 
        return p**2*pke*(nZ*sig)*epf/(p*qval)*2*np.pi
    return 0

if fg != 1:
    with open('pke/pke12_tot.data') as f:
        n_e, n_p = f.readline().split()
        n_e = int(n_e)
        n_p = int(n_p)
        p = np.empty(n_p)
        pke = np.empty((n_e,n_p))
        dp = np.empty(n_p)
        xe = np.empty(n_e)
        for j in range(n_p):
            p[j] = float(f.readline())
            for i in range(int(n_e/4)):
                tokens = f.readline().split()
                for k in range(4):
                    xe[4*i+k] = tokens[2*k]
                    pke[4*i+k,j] = tokens[2*k+1]

pmax = p[-1]
hp = p[1]-p[0]
he = xe[1]-xe[0]

norm = 0.0
for i in range(n_p):
    dp[i] = sum(pke[:,i])*he
norm = sum(p**2*dp)*4*np.pi*hp
pke /= norm
print('n(k) norm initial = {0}'.format(norm))

pke_interp = interpolate.interp2d(p,xe,pke,kind='linear')

omega = []
mom = []
energies = []
p_px = []
p_py = []
p_pz = []
p_E = []
wgts = []
escape = []

def GenerateEvent(x):
    dw = wmax - wmin
    dp = p[-1] - p[0]
    de = xe[-1] - xe[0]
    w = dw*x[0] + wmin
    p_int = dp*x[1] + p[0]
    e_int = de*x[2] + xe[0]
    eef = ee - w
    Q2 = 2.0*ee*eef*(1.0-coste)

    qval = np.sqrt(Q2 + w**2)
    f_o = f_eval(nZ,p_int,e_int,w,qval,mqe,pke_interp(p_int,e_int),fg,kf,thetalept,ee,iform)

    wgt = f_o*1e9*dw*dp*de

    return wgt

integ = vegas.Integrator([[0,1],[0,1],[0,1]])

result = integ(GenerateEvent, nitn=10, neval=1000)

print(result.summary())

from nuChic.FourVector import Vec4
from nuChic.Nucleus import Nucleus
from nuChic.Constants import MeV, GeV, fm, mN
from nuChic.Cascade import FSI

argon_nucleus = Nucleus(6,12, 8.6*MeV, 225*MeV)
fsi = FSI(argon_nucleus, 1)


nitn = 2
count = 0
for i in range(nitn):
    for x, weight in integ.random():
    
        if count % 1000 == 0:
            print(count)
        count += 1

        dw = wmax - wmin
        dp = p[-1] - p[0]
        de = xe[-1] - xe[0]
        w = dw*x[0] + wmin
        p_int = dp*x[1] + p[0]
        e_int = de*x[2] + xe[0]
        eef = ee - w
        Q2 = 2.0*ee*eef*(1.0-coste)

        qval = np.sqrt(Q2 + w**2)
        f_o = f_eval(nZ,p_int,e_int,w,qval,mqe,pke_interp(p_int,e_int),fg,kf,thetalept,ee,iform)
   
        if f_o == 0:
            continue

        wgt = f_o*1e9*dw*dp*de
   
        ep = np.sqrt(p_int**2+mqe**2)
        wt = w - e_int + mqe - ep
        cost_te = ((wt+ep)**2 - p_int**2 - qval**2 - mqe**2)/(2.0*p_int*qval)
        pf = np.sqrt(p_int**2 + qval**2 + 2*qval*p_int*cost_te)
        epf = np.sqrt(mqe**2+pf**2)
        phi = 2.0*np.pi*np.random.random()

        xq = qval/hbarc
        xk = p_int/hbarc
        xp = pf/hbarc

        q2 = xq**2
        p2 = xk**2
        pf2 = xp**2
        cosa = ((pf2-p2-q2)/2.0/xk/xq)
        sina2 = 1-cosa**2

        omega.append(w)
        mom.append(p_int)
        energies.append(e_int)
        p_E.append(epf)
        p_px.append(xk*np.sqrt(sina2)*np.cos(phi)*hbarc)
        p_py.append(xk*np.sqrt(sina2)*np.sin(phi)*hbarc)
        p_pz.append(xk*cosa*hbarc + qval)
        wgts.append(wgt[0]*weight/nitn)

        momentum = Vec4(epf*MeV, pf*np.sqrt(sina2)*np.cos(phi)*MeV, pf*np.sqrt(sina2)*np.sin(phi)*MeV, pf*cosa*MeV)
        fsi.kick(momentum)
        escape.append(len(fsi()))
        fsi.reset()

print(sum(wgts))
print(count)

import pandas as pd

data = pd.read_csv('../data/12C.dat',header=None,sep='\s+',names=['Z','A','Ee','Angle','omega','dsigma','error','citation'])

energy = 0.730
angle = 37.1

maskEnergy = data['Ee'] == energy
maskAngle = data['Angle'] == angle
data_tmp = data[maskEnergy & maskAngle]
data_omega = data_tmp['omega'].values
data_dsigma = data_tmp['dsigma'].values
data_error = data_tmp['error'].values

# Load theory prediction
filename = 'C12_{0}_{1}p{2}.out'.format(int(energy*1000),int(angle),int((angle*10)%10))
theory = pd.read_csv(filename,header=None,names=['omega','dsigma','error'],sep='\s+')
omega_f = theory['omega'].values
dsigma = theory['dsigma'].values

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1)
ax1.errorbar(data_omega*1000, data_dsigma, yerr=data_error, fmt='o', color='red')
ax1.plot(omega_f,dsigma, color='red', ls='steps')
ax1.hist(omega,weights=wgts,bins=645,color='blue')
ax2.hist(mom,weights=wgts,bins=400)
ax3.hist(energies,weights=wgts,bins=400)

fig2, ((ax4, ax5), (ax6, ax7)) = plt.subplots(nrows=2, ncols=2)
ax4.hist(p_E,weights=wgts,bins=400)
ax5.hist(p_px,weights=wgts,bins=400)
ax6.hist(p_py,weights=wgts,bins=400)
ax7.hist(p_pz,weights=wgts,bins=400)

fig3, ax7 = plt.subplots(nrows=1, ncols=1)
ax7.hist(escape,weights=wgts,bins=np.linspace(0,10,11))

plt.show()


