""" Implement inclusive cross-section calculations"""
import os

import xsec
import numpy as np
from scipy import interpolate
from absl import flags
from absl import logging

from .four_vector import Vec4
from .nucleus import Nucleus
from .constants import MEV, HBARC, MQE
from .folding import Folding

FLAGS = flags.FLAGS
flags.DEFINE_bool(
    'pwia', None, 'Flag to turn on/off the plane-wave impluse approximation')

DIR, FILE = os.path.split(__file__)


class Inclusive:
    """ Base class to implement inclusive cross-section calculations."""
    def __init__(self, nucleus, energy, thetalept):
        self.beam_energy = energy
        self.thetalept = thetalept*np.pi/180.0

        if not isinstance(nucleus, Nucleus):
            raise ValueError('Requires Nucleus as input, got {}'.format(
                type(nucleus)))

        self.nucleus = nucleus

        self.wrange = [1, 650]

        self.folding = None
        if FLAGS.folding:
            self.folding = Folding()

    @property
    def n_z(self):
        """ Get the number of protons in the nucleus. """
        return self.nucleus.protons

    @property
    def k_f(self):
        """ Get the Fermi momentum. """
        return self.nucleus.kf

    @property
    def wmin(self):
        """ Get the minimum allowed omega. """
        return self.wrange[0]

    @property
    def wmax(self):
        """ Get the maximum allowed omega. """
        return self.wrange[1]

    @property
    def coste(self):
        """ Get the cosine of the outgoing lepton angle. """
        return np.cos(self.thetalept)

    def generate_weight(self, point):
        """ Generate weight for a given point. """
        raise NotImplementedError()

    def generate_momentum(self, point):
        """ Generate momentum for a given point. """
        raise NotImplementedError()


class Quasielastic(Inclusive):
    """ Class to calculate quasielastic scattering of an electron
    on a nucleus."""
    def __init__(self, nucleus, energy, thetalept, fg=0):
        logging.info('Initializing Quasielastic calculation')

        super().__init__(nucleus, energy, thetalept)

        self.f_g = fg
        self.iform = 2
        xsec.dirac_matrices.dirac_matrices_in(MQE/HBARC)

        if fg != 1:
            with open(os.path.join(DIR, 'pke', 'pke12_tot.data')) as infile:
                n_e, n_p = infile.readline().split()
                n_e = int(n_e)
                n_p = int(n_p)
                self.mom = np.empty(n_p)
                self.pke = np.empty((n_e, n_p))
                d_p = np.empty(n_p)
                self.energy = np.empty(n_e)
                for j in range(n_p):
                    self.mom[j] = float(infile.readline())
                    for i in range(int(n_e/4)):
                        tokens = infile.readline().split()
                        for k in range(4):
                            self.energy[4*i+k] = tokens[2*k]
                            self.pke[4*i+k, j] = tokens[2*k+1]

        d_p = np.sum(self.pke, axis=0)*(self.energy[1]-self.energy[0])
        norm = np.sum(self.mom**2*d_p, axis=-1)*4*np.pi*(
            self.mom[1]-self.mom[0])
        self.pke /= norm
        logging.info('n(k) norm initial = {0}'.format(norm))

        self.pke = interpolate.interp2d(self.mom, self.energy,
                                        self.pke, kind='linear')

    @property
    def pmin(self):
        """ Minimum allowed momentum in initial state. """
        return self.mom[0]

    @property
    def pmax(self):
        """ Maximum allowed momentum in initial state. """
        return self.mom[-1]

    @property
    def emin(self):
        """ Minimum allowed energy in initial state. """
        return self.energy[0]

    @property
    def emax(self):
        """ Maximum allowed energy in initial state. """
        return self.energy[-1]

    def _eval(self, mom, energy, omega, qval, pke):
        e_out = np.sqrt(mom**2+MQE**2)
        if self.f_g == 1:
            omegat = omega
        else:
            omegat = omega-energy+MQE-e_out

        u_pq = 0.0
        if FLAGS.folding:
            tkin_pf = np.sqrt(qval**2+MQE**2)-MQE
            if self.folding.kin_max < tkin_pf < self.folding.kin_min:
                u_pq = self.folding.kinematic(tkin_pf)

        cost_te = (((omegat+e_out+u_pq)**2 - mom**2 - qval**2 - MQE**2)
                   / (2*mom*qval))
        if abs(cost_te) > 1:
            logging.debug('omegat = %e, ep = %e, p = %e, qval = %e, mqe = %e'
                          % (omegat, e_out, mom, qval, MQE))
            return 0

        mom_f = np.sqrt(mom**2 + qval**2 + 2*qval*mom*cost_te)
        epf = np.sqrt(MQE**2 + mom_f**2)
        if mom_f >= self.k_f:
            phi = 2.0*np.pi*np.random.rand(1)
            sig = xsec.cc1(qval/HBARC, omega, omegat, mom/HBARC, mom_f/HBARC,
                           phi, self.energy, self.thetalept, self.iform)
            return mom**2*pke[0]*(self.n_z*sig)*epf/(mom*qval)*2*np.pi
        return 0

    def generate_weight(self, point):
        domega = self.wmax - self.wmin
        omega = domega*point[0] + self.wmin
        denergy = self.emax - self.emin
        e_int = denergy*point[1] + self.emin

        dmom = self.pmax - self.pmin
        p_int = dmom*point[2] + self.pmin
        if FLAGS.folding:
            domegap = self.wmax - self.wmin
            omegap = domegap*point[3] + self.wmin

        eef = self.beam_energy - omega
        q_2 = 2.0 * self.beam_energy * eef * (1.0 - self.coste)
        qval = np.sqrt(q_2 + omega**2)

#        dcos = np.sqrt(mqe**2+self.qval**2-(self.e_int+self.w)**2)/self.qval
#        print(dcos)
#        cost_te = 2*dcos*x[2]-dcos

        wgt = self._eval(p_int, e_int, omega, qval, self.pke(p_int, e_int))
        wgt *= 1e9*domega*dmom*denergy

        if FLAGS.folding:
            qval_f = np.sqrt(2.0*self.energy*(self.energy-omegap)
                             * (1.0-self.coste) + omegap**2)

            wgt_f = self._eval(p_int, e_int, omegap, qval_f,
                               self.pke(p_int, e_int))
            wgt_f *= 1e9*domegap*dmom*denergy
            wgt_f = self.folding(omega, omegap)*wgt_f * \
                domega + self.folding.transparency*wgt

            return wgt_f, wgt

        return wgt

    def generate_momentum(self, point):
        domega = self.wmax - self.wmin
        omega = domega*point[0] + self.wmin
        denergy = self.emax - self.emin
        e_int = denergy*point[1] + self.emin

        dmom = self.pmax - self.pmin
        p_int = dmom*point[2] + self.pmin

        eef = self.beam_energy - omega
        q_2 = 2.0 * self.beam_energy * eef * (1.0 - self.coste)
        qval = np.sqrt(q_2 + omega**2)

        e_out = np.sqrt(p_int**2+MQE**2)
        omegat = omega - e_int + MQE - e_out
        cost_te = ((omegat+e_out)**2 - p_int**2 - qval **
                   2 - MQE**2)/(2.0*p_int*qval)
        if abs(cost_te) > 1:
            return None
        mom_f = np.sqrt(p_int**2 + qval**2 +
                        2*qval*p_int*cost_te)
        epf = np.sqrt(MQE**2+mom_f**2)
        phi = 2.0*np.pi*np.random.random()

        x_q = qval/HBARC
        x_k = p_int/HBARC
        x_p = mom_f/HBARC

        q_2 = x_q**2
        p_2 = x_k**2
        pf2 = x_p**2
        cosa = ((pf2-p_2-q_2)/2.0/x_k/x_q)
        sina2 = 1-cosa**2

        momentum = Vec4(epf*MEV, mom_f*np.sqrt(sina2)*np.cos(phi)*MEV,
                        mom_f*np.sqrt(sina2)*np.sin(phi)*MEV, mom_f*cosa*MEV)

        return momentum
