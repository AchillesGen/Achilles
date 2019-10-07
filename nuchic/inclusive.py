""" Implement inclusive cross-section calculations"""

import xsec
import numpy as np
from scipy import interpolate
from absl import flags
from absl import logging

from .four_vector import Vec4
from .nucleus import Nucleus
from .constants import MEV, HBARC, MQE
from .folding import Folding
from .config import SETTINGS
from .utils import make_path

FLAGS = flags.FLAGS
flags.DEFINE_bool(
    'pwia', None, 'Flag to turn on/off the plane-wave impluse approximation')


class Inclusive:
    """ Base class to implement inclusive cross-section calculations."""
    def __init__(self, nucleus, energy, thetalept):
        self.beam_energy = energy
        self.thetalept = thetalept*np.pi/180.0

        if not isinstance(nucleus, Nucleus):
            raise ValueError('Requires Nucleus as input, got {}'.format(
                type(nucleus)))

        self.nucleus = nucleus

        self.wrange = [1, energy]

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

    def generate_momentum(self, variables, qval):
        """ Generate momentum for a given point. """
        raise NotImplementedError()


class Quasielastic(Inclusive):
    """ Class to calculate quasielastic scattering of an electron
    on a nucleus."""
    def __init__(self, fg, *args, **kwargs):
        logging.info('Initializing Quasielastic calculation')

        super().__init__(*args, **kwargs)

        self.width = 1e3
        self.f_g = fg
        self.iform = 2
        xsec.dirac_matrices.dirac_matrices_in(MQE/HBARC)

        if fg != 1:
            with open(make_path('pke12_tot.data', 'pke')) as infile:
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

    def _eval(self, variables, qval, pke):
        e_out = np.sqrt(variables['mom']**2+MQE**2)
        if self.f_g == 1:
            omegat = variables['omega']
        else:
            omegat = variables['omega']-variables['energy']+MQE-e_out

        u_pq = 0.0
        if FLAGS.folding:
            tkin_pf = np.sqrt(qval**2+MQE**2)-MQE
            if self.folding.kin_max < tkin_pf < self.folding.kin_min:
                u_pq = self.folding.kinematic(tkin_pf)

        cost_te = (((omegat+e_out+u_pq)**2 - variables['mom']**2 - qval**2
                    - MQE**2)
                   / (2*variables['mom']*qval))
        if abs(cost_te) > 1:
            logging.debug('omegat = %e, ep = %e, p = %e, qval = %e, mqe = %e'
                          % (omegat, e_out, variables['mom'], qval, MQE))
            return 0
        variables['cost_te'] = cost_te

        mom_f = np.sqrt(variables['mom']**2 + qval**2
                        + 2*qval*variables['mom']*variables['cost_te'])
        if mom_f >= self.k_f:
            phi = 2.0*np.pi*np.random.rand(1)
            sig = xsec.cc1(qval/HBARC, variables['omega'], omegat,
                           variables['mom']/HBARC,
                           mom_f/HBARC, phi, self.beam_energy, self.thetalept,
                           self.iform)
            return (variables['mom']**2*pke[0]*(self.n_z*sig)
                    * np.sqrt(MQE**2 + mom_f**2)*2*np.pi
                    / (variables['mom']*qval))
        return 0

    # TODO: How to best implement this numerically?
    # Currently, this is very unstable from run to run
    def _delta(self, omegap, mom, qval, cost):
        arg = omegap**2 - mom**2 - qval**2 - MQE**2 - 2*mom*qval*cost
        return 1.0/(self.width * np.sqrt(np.pi))*np.exp(-(arg/self.width)**2)

    def _map_vars(self, point):
        domega = self.wmax - self.wmin
        omega = domega*point[0] + self.wmin
        denergy = self.emax - self.emin
        e_int = denergy*point[1] + self.emin
        dmom = self.pmax - self.pmin
        p_int = dmom*point[2] + self.pmin

        variables = {'mom': p_int,
                     'energy': e_int,
                     'omega': omega}

        ps_wgt = domega*denergy*dmom
        qval = np.sqrt((2.0 * self.beam_energy * (self.beam_energy - omega)
                        * (1.0 - self.coste)) + omega**2)

        if FLAGS.folding:
            domegap = self.wmax - self.wmin
            omegap = domegap*point[3] + self.wmin
            variables['omegap'] = omegap
            ps_wgt = domegap*denergy*dmom
            return variables, ps_wgt, qval, domega

        return variables, ps_wgt, qval, None

    def generate_weight(self, point):
        variables, ps_wgt, qval, domega = self._map_vars(point)

        wgt = self._eval(variables, qval, self.pke(variables['mom'],
                                                   variables['energy']))
        wgt *= 1e9*ps_wgt

        if FLAGS.folding:
            qval_f = np.sqrt(2.0*self.beam_energy*(self.beam_energy
                                                   - variables['omegap'])
                             * (1.0-self.coste) + variables['omegap']**2)

            variables['omega'], variables['omegap'] = (variables['omegap'],
                                                       variables['omega'])

            wgt_f = self._eval(variables, qval_f,
                               self.pke(variables['mom'], variables['energy']))
            wgt_f *= 1e9*ps_wgt
            wgt_f = self.folding(SETTINGS.run['folding_func'],
                                 variables['omega'],
                                 variables['omegap'])\
                * wgt_f * domega + self.folding.transparency*wgt

            variables['omega'], variables['omegap'] = (variables['omegap'],
                                                       variables['omega'])

            return wgt_f, variables, qval

        return wgt, variables, qval

    def generate_momentum(self, variables, qval):
        e_out = np.sqrt(variables['mom']**2+MQE**2)
        if self.f_g == 1:
            omegat = variables['omega']
        else:
            omegat = variables['omega']-variables['energy']+MQE-e_out
        u_pq = 0.0
        if FLAGS.folding:
            tkin_pf = np.sqrt(qval**2+MQE**2)-MQE
            if self.folding.kin_max < tkin_pf < self.folding.kin_min:
                u_pq = self.folding.kinematic(tkin_pf)

        cost_te = (((omegat+e_out+u_pq)**2 - variables['mom']**2 - qval**2
                    - MQE**2)
                   / (2*variables['mom']*qval))
        if abs(cost_te) > 1:
            logging.debug('omegat = %e, ep = %e, p = %e, qval = %e, mqe = %e'
                          % (omegat, e_out, variables['mom'], qval, MQE))
            return None
        variables['cost_te'] = cost_te
        mom_f = np.sqrt(variables['mom']**2 + qval**2 +
                        2*qval*variables['mom']*variables['cost_te'])
        epf = np.sqrt(MQE**2+mom_f**2)
        phi = 2.0*np.pi*np.random.random()

        x_q = qval/HBARC
        x_k = variables['mom']/HBARC
        x_p = mom_f/HBARC

        q_2 = x_q**2
        p_2 = x_k**2
        pf2 = x_p**2
        cosa = ((pf2-p_2-q_2)/2.0/x_k/x_q)

        momentum = Vec4(epf*MEV,
                        mom_f*np.sqrt(1-cosa**2)*np.cos(phi)*MEV,
                        mom_f*np.sqrt(1-cosa**2)*np.sin(phi)*MEV,
                        mom_f*cosa*MEV)

        return momentum
