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
ee = 1511.0
thetalept = 89.97
nZ = 1
nA = 2
fg = 0
kf = 225
nw = 400
wmin = 800
wmax = 1200

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

# Create new interpolation function for vegas to use after normalization
pke_interp = interpolate.interp1d(p0,pke0[0],kind='cubic')

omega = []
mom = []
p_px = []
p_py = []
p_pz = []
p_E = []
wgts = []
escape = []

def GenerateEvent(x):
    dw = wmax-wmin
    dp = p[-1] - p[0]
    w = dw*x[0]+wmin
    p_int = dp*x[1]+p[0]
    eef = ee - w
    Q2 = 2.0*ee*eef*(1.0 - coste)

    # Compute the response functions in the impulse approximation
    qval = np.sqrt(Q2 + w**2)
    f_o = f_eval(nZ,p_int,w,qval,mqe,pke_interp(p_int),fg,kf,thetalept,ee,iform)

    wgt = f_o*1e9*dw*dp

#    if fill:
#        omega.append(w)
#        mom.append(p_int)
#        wgts.append(wgt/1e9)

    return wgt

integ = vegas.Integrator([[0,1],[0,1]])

#fill = False
integ(GenerateEvent,nitn=10,neval=1e5)

#fill = True
#result = integ(GenerateEvent,nitn=10,neval=1e6)
# Print a summary of the VEGAS Integration results
#print("Summary of VEGAS integration results")
#print(result.summary())
#print("The best result coming from VEGAS is {0} pb".format(result))
#print(sum(wgts))

from nuChic.FourVector import Vec4
import numpy as np
from nuChic.Nucleus import Nucleus
from nuChic.Constants import hbarc, MeV, GeV, fm, mN
from nuChic.Cascade import FSI

argon_nucleus = Nucleus(6,12, 8.6*MeV, 225*MeV)
fsi = FSI(argon_nucleus, 1)

nitn = 1
for i in range(nitn):
    for x, weight in integ.random():
        dw = wmax-wmin
        dp = p[-1] - p[0]
        w = dw*x[0]+wmin
        p_int = dp*x[1]+p[0]
        eef = ee - w
        Q2 = 2.0*ee*eef*(1.0 - coste)
    
        # Compute the response functions in the impulse approximation
        qval = np.sqrt(Q2 + w**2)
        f_o = f_eval(nZ,p_int,w,qval,mqe,pke_interp(p_int),fg,kf,thetalept,ee,iform)
   
        if f_o == 0:
            continue

        wgt = f_o*1e9*dw*dp
   
        ep = np.sqrt(p_int**2+mqe**2)
        wt = w - 2.0
        cost_te = ((wt+ep)**2 - p_int**2 - qval**2 - mqe**2)/(2.0*p_int*qval)
        pf = np.sqrt(p_int**2 + qval**2 + 2*qval*p_int*cost_te)
        epf = np.sqrt(mqe**2+pf**2)
        phi = 2.0*pi*np.random.random()

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
        p_E.append(epf)
        p_px.append(xk*np.sqrt(sina2)*np.cos(phi))
        p_py.append(xk*np.sqrt(sina2)*np.sin(phi))
        p_pz.append(xk*cosa)
        wgts.append(wgt*weight/nitn)

        momentum = Vec4(epf, xk*np.sqrt(sina2)*np.cos(phi), xk*np.sqrt(sina2)*np.sin(phi), xk*cosa)
        fsi.kick(momentum)
        escape.append(len(fsi()))
        fsi.reset()


print(sum(wgts))

import pandas as pd

data = pd.read_csv('../../FNALNeuGen/data/2H.dat',header=None,sep='\s+',names=['Z','A','Ee','Angle','omega','dsigma','error','citation'])

energy = 1.511
angle = 89.97

maskEnergy = data['Ee'] == energy
maskAngle = data['Angle'] == angle
data_tmp = data[maskEnergy & maskAngle]
data_omega = data_tmp['omega'].values
data_dsigma = data_tmp['dsigma'].values
data_error = data_tmp['error'].values

# Load theory prediction
filename = '../../code_deuteron/H2_{0}_{1}p{2}_IIb.out'.format(int(energy*1000),int(angle),int((angle*10)%10))
theory = pd.read_csv(filename,header=None,names=['omega','dsigma','error'],sep='\s+')
omega_f = theory['omega'].values
dsigma = theory['dsigma'].values

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
ax1.errorbar(data_omega*1000, data_dsigma, yerr=data_error, fmt='.')
ax1.plot(omega_f,dsigma, color='red', ls='steps')
ax1.hist(omega,weights=wgts,bins=400)
ax2.hist(mom,weights=wgts,bins=400)

fig2, ((ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=2, ncols=2)
ax3.hist(p_E,weights=wgts,bins=400)
ax4.hist(p_px,weights=wgts,bins=400)
ax5.hist(p_py,weights=wgts,bins=400)
ax6.hist(p_pz,weights=wgts,bins=400)

fig3, ax7 = plt.subplots(nrows=1, ncols=1)
ax7.hist(escape,weights=wgts,bins=400)

plt.show()


