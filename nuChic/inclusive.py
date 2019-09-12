import xsec
import numpy as np
from scipy import interpolate
from absl import flags
from absl import logging
import os

from .FourVector import Vec4
from .Nucleus import Nucleus
from .Constants import MeV, GeV, fm, hbarc, mqe
from .Cascade import FSI
from .Folding import Folding

FLAGS = flags.FLAGS
flags.DEFINE_bool('pwia',None,'Flag to turn on/off the plane-wave impluse approximation')

DIR, FILE = os.path.split(__file__)

class Inclusive:
    def __init__(self):
        pass

class Quasielastic(Inclusive):
    def __init__(self,nucleus,ee,thetalept,fg=0):
        logging.info('Initializing Quasielastic calculation')
        
        self.ee = ee
        self.thetalept = thetalept
        self.nZ = nucleus.Z
        self.kf = nucleus.kf/MeV
        self.fg = fg
        self.wmin = 1
        self.wmax = 650
        self.iform = 2
        self.thetalept*=np.pi/180.0
        self.coste = np.cos(thetalept)
        xsec.dirac_matrices.dirac_matrices_in(mqe/hbarc)

        self.folding = None
        if FLAGS.folding:
            self.folding = Folding()

        if fg != 1:
            with open(os.path.join(DIR,'pke','pke12_tot.data')) as f:
                n_e, n_p = f.readline().split()
                n_e = int(n_e)
                n_p = int(n_p)
                self.p = np.empty(n_p)
                pke = np.empty((n_e,n_p))
                dp = np.empty(n_p)
                self.xe = np.empty(n_e)
                for j in range(n_p):
                    self.p[j] = float(f.readline())
                    for i in range(int(n_e/4)):
                        tokens = f.readline().split()
                        for k in range(4):
                            self.xe[4*i+k] = tokens[2*k]
                            pke[4*i+k,j] = tokens[2*k+1]

        hp = self.p[1] - self.p[0]
        he = self.xe[1] - self.xe[0]

        dp = np.sum(pke,axis=0)*he
        norm = np.sum(self.p**2*dp,axis=-1)*4*np.pi*hp
        pke /= norm
        logging.info('n(k) norm initial = {0}'.format(norm))

        self.pke = interpolate.interp2d(self.p, self.xe, pke, kind='linear')

    @property
    def pmin(self):
        return self.p[0]

    @property
    def pmax(self):
        return self.p[-1]

    @property
    def emin(self):
        return self.xe[0]

    @property
    def emax(self):
        return self.xe[-1]

    def _eval(self,p,e,w,qval,pke):
        ep = np.sqrt(p**2+mqe**2)
        if self.fg == 1:
            wt = w
        else:
            wt = w-e+mqe-ep

        u_pq = 0.0
        #if FLAGS.folding:
        #    tkin_pf = np.sqrt(qval**2+mqe**2)-mqe
        #    if tkin_pf < self.folding.kin_min and tkin_pf > self.folding.kin_max:
        #        u_pq = self.folding.kinematic(tkin_pf)

        cost_te = ((wt+ep+u_pq)**2 - p**2 - qval**2 - mqe**2)/(2*p*qval)
        if abs(cost_te) > 1:
            logging.debug('wt = %e, ep = %e, p = %e, qval = %e, mqe = %e' 
                % (wt,ep,p,qval,mqe))
            return 0

        pf = np.sqrt(p**2 + qval**2 + 2*qval*p*cost_te)
        epf = np.sqrt(mqe**2 + pf**2)
        if pf >= self.kf:
            phi = 2.0*np.pi*np.random.rand(1)
            sig = xsec.cc1(qval/hbarc,w,wt,p/hbarc,pf/hbarc,phi,self.ee,
                self.thetalept,self.iform) 
            return p**2*pke[0]*(self.nZ*sig)*epf/(p*qval)*2*np.pi
        return 0

    def GenerateWeight(self,x):
        dw = self.wmax - self.wmin
        self.w = dw*x[0] + self.wmin
#        mass = self.ee**2*(1-self.coste)/(mqe + self.ee*(1-self.coste))
#        mass2 = mass**2
#        ymax = np.arctan((0-mass2)/mass2)
#        ymin = np.arctan((self.ee**2-mass2)/mass2)
#        self.w = np.sqrt(mass2+mass2*np.tan(ymin+x[0]*(ymax-ymin)))
#        dw = (ymin - ymax)/(mass2/((self.w**2-mass2)**2+mass2**2))/(2*self.w)
        de = self.w - self.emin
        self.e_int = de*x[1] + self.emin

        dp = self.pmax - self.pmin
        self.p_int = dp*x[2] + self.pmin
        if FLAGS.folding:
#            mass = self.ee**2*(1-self.coste)/(mqe + self.ee*(1-self.coste))
#            mass2 = mass**2
#            ymax = np.arctan((0-mass2)/mass2)
#            ymin = np.arctan((self.ee**2-mass2)/mass2)
#            self.wp = np.sqrt(mass2+mass2*np.tan(ymin+x[3]*(ymax-ymin)))
#            dwp = (ymin - ymax)/(mass2/((self.wp**2-mass2)**2+mass2**2))/(2*self.wp)
            dwp = self.wmax - self.wmin
            self.wp = dwp*x[3] + self.wmin

        eef = self.ee - self.w
        Q2 = 2.0*self.ee*eef*(1.0-self.coste)

        self.qval = np.sqrt(Q2 + self.w**2)

#        dcos = np.sqrt(mqe**2+self.qval**2-(self.e_int+self.w)**2)/self.qval
#        print(dcos)
#        cost_te = 2*dcos*x[2]-dcos
#
#        print(de,self.w,self.e_int,cost_te,-mqe**2+(-1+cost_te**2)*self.qval**2 + (mqe-self.e_int+self.w)**2)
#        self.p_int = -cost_te * self.qval - np.sqrt(-mqe**2+(-1+cost_te**2)*self.qval**2 + (mqe-self.e_int+self.w)**2) 
#
#        print('minus: ', self.p_int)
#
#        self.p_int = -cost_te * self.qval + np.sqrt(-mqe**2+(-1+cost_te**2)*self.qval**2 + (mqe-self.e_int+self.w)**2) 
#
#        print('plus: ', self.p_int)

        pke = self.pke(self.p_int, self.e_int)
        wgt = self._eval(self.p_int,self.e_int,self.w,self.qval,pke)

        wgt *= 1e9*dw*dp*de

        if FLAGS.folding:
            eef_f = self.ee - self.wp
            Q2_f = 2.0*self.ee*eef_f*(1.0-self.coste)
            qval_f = np.sqrt(Q2_f + self.wp**2)

            wgt_f = self._eval(self.p_int,self.e_int,self.wp,qval_f,pke)
            wgt_f *= 1e9*dwp*dp*de
            wgt_f = self.folding(self.w,self.wp)*wgt_f*dw + self.folding.TA*wgt

            return wgt_f, wgt

        return wgt

    def GenerateMomentum(self):
        ep = np.sqrt(self.p_int**2+mqe**2)
        wt = self.w - self.e_int + mqe - ep
        cost_te = ((wt+ep)**2 - self.p_int**2 - self.qval**2 - mqe**2)/(2.0*self.p_int*self.qval)
        if abs(cost_te) > 1:
            return None
        pf = np.sqrt(self.p_int**2 + self.qval**2 + 2*self.qval*self.p_int*cost_te)
        epf = np.sqrt(mqe**2+pf**2)
        phi = 2.0*np.pi*np.random.random()
    
        xq = self.qval/hbarc
        xk = self.p_int/hbarc
        xp = pf/hbarc
    
        q2 = xq**2
        p2 = xk**2
        pf2 = xp**2
        cosa = ((pf2-p2-q2)/2.0/xk/xq)
        sina2 = 1-cosa**2
        
        momentum = Vec4(epf*MeV, pf*np.sqrt(sina2)*np.cos(phi)*MeV, 
            pf*np.sqrt(sina2)*np.sin(phi)*MeV, pf*cosa*MeV)

        return momentum
