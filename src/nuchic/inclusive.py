""" Implement inclusive cross-section calculations"""

from collections import namedtuple
from scipy import interpolate
# from absl import flags
# from absl import logging
import numpy as np

# import xsec
from .four_vector import Vec4
from .constants import MEV, HBARC, MQE, TO_NB
from .folding import Folding
from .config import settings
from .utils import make_path

# FLAGS = flags.FLAGS
# flags.DEFINE_bool(
#     'pwia', None, 'Flag to turn on/off the plane-wave impluse approximation')


class Inclusive:
    """ Base class to implement inclusive cross-section calculations."""
    def __init__(self, energy, thetalept):
        self.beam_energy = energy
        self.thetalept = thetalept*np.pi/180.0

        self.nucleus = settings().nucleus

        self.wrange = [1, energy]

        self.folding = None
        # if FLAGS.folding:
        #     self.folding = Folding()

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

    def generate_momentum(self, var):
        """ Generate momentum for a given point. """
        raise NotImplementedError()


def read_spectral_function_data():
    """
    Reads spectral function data from the file 'pke12_tot.data'.
    The data in this file appears in several stanzas.
    The first line of the file specifies the sizes of the stanzas

    has the format:
    n_e, n_p
    energy (or is it momentum??)
    energy, the value of the spectral function


    """
    spectral_params = namedtuple(
        'SpectralParams',
        ['n_e', 'n_p', 'mom', 'energy', 'pke', 'd_p', 'norm']
    )
    with open(make_path('pke12_tot.data', 'pke')) as infile:
        n_e, n_p = infile.readline().split()
        n_e = int(n_e)
        n_p = int(n_p)
        mom = np.empty(n_p)
        pke = np.empty((n_e, n_p))
        d_p = np.empty(n_p)
        energy = np.empty(n_e)
        for j in range(n_p):
            mom[j] = float(infile.readline())
            for i in range(int(n_e/4)):
                tokens = infile.readline().split()
                for k in range(4):
                    energy[4*i+k] = tokens[2*k]
                    pke[4*i+k, j] = tokens[2*k+1]
    d_p = np.sum(pke, axis=0) * (energy[1] - energy[0])
    norm = np.sum(mom**2 * d_p, axis=-1) * 4 * np.pi * (mom[1] - mom[0])
    pke /= norm
    return spectral_params(n_e, n_p, mom, energy, pke, d_p, norm)


Variables = namedtuple(
    'Variables',
    ['mom', 'energy', 'omega', 'omegap', 'cost_te', 'qval', 'qval_f']  #,
    # defaults=[None, None, None, None, None, None, None]
)


def _update_variables(var, attr, new_val):
    """
    Creates a new instance of the named tuple Variable with an updated value
    for the specified attributed.
    Args:
        var: Variables
        attr: str, the name of the attribute to update
        new_val: the value of the attribute for updating
    Returns:
        Variables
    """
    attrs = var._fields
    kwargs = {attr_i: getattr(var, attr_i) for attr_i in attrs}
    kwargs[attr] = new_val
    return Variables(**kwargs)


class Quasielastic(Inclusive):
    """
    Class to calculate quasielastic scattering of an electron on a nucleus.
    Args:
        fg: int, flag for whether or not use the spectral function
    Raises:
        NotImplementedError, when fg == 1.
    """
    def __init__(self, fg, *args, **kwargs):
        logging.info('Initializing Quasielastic calculation')

        super(Quasielastic, self).__init__(*args, **kwargs)
        self.f_g = fg
        self.width = 1e3  # Width of delta function, [width] = MeV^2
        self.iform = 2  # hard-coded flag which gets passed to Noemi's code
        xsec.dirac_matrices.dirac_matrices_in(MQE/HBARC)

        if fg == 1:
            raise NotImplementedError(
                'Quasielastic only implemented for fg != 1.')
        spectral_params = read_spectral_function_data()

        self.mom = spectral_params.mom
        self.energy = spectral_params.energy
        self.pke = spectral_params.pke
        self.d_p = spectral_params.d_p

        logging.info('n(k) norm initial = {0}'.format(spectral_params.norm))

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

    def _eval(self, var, pke):
        """
        Evaluates the weight using the charced-current cross secion and the
        nuclear physics constraints on the process.
        Args:
            var: namedtuple Variables containing kinematic quantities
            pke: the spectral function

        Returns:
            (wgt, var): the wgt and the kinematic variables, which are updated
                to include 'cost_te' from the energy-conserving delta function.
        """
        # Compute outgoing energy from relativistic dispersion
        e_out = np.sqrt(var.mom**2 + MQE**2)
        if self.f_g == 1:
            omegat = var.omega
        else:
            omegat = var.omega - var.energy + MQE - e_out

        # Adjust the momentum due to the optical potential
        u_pq = 0.0
        # if FLAGS.folding:
        #     tkin_pf = np.sqrt(var.qval**2 + MQE**2) - MQE
        #     if self.folding.kin_min < tkin_pf < self.folding.kin_max:
        #         u_pq = self.folding.kinematic(tkin_pf)

        # Enforce value of cos(theta) from the energy-conserving delta function
        cost_te = (((omegat + e_out - u_pq)**2 - var.mom**2 - var.qval**2 - MQE**2)
                   / (2 * var.mom * var.qval))
        if abs(cost_te) > 1:
            # Something bad happened, why is |cos(theta)| > 1?
            logging.debug('omegat = %e, ep = %e, p = %e, qval = %e, mqe = %e'
                          % (omegat, e_out, var.mom, var.qval, MQE))
            return np.array([0.0, 0.0]), var  # vanishing wgt

        var = _update_variables(var, 'cost_te', cost_te)

        # Compute the final momentum
        mom_f = np.sqrt(var.mom**2 + var.qval**2
                        + 2 * var.qval * var.mom * var.cost_te)

        # Check if the final momentum falls below the Fermi momentum
        if mom_f < self.k_f:
            return np.array([0.0, 0.0]), var  # vanishing wgt

        phi = 2.0 * np.pi * np.random.rand(1)  # Random azimuthal angle
        # Compute charged-current cross section
        sig = xsec.cc1(var.qval/HBARC, var.omega, omegat,
                       var.mom/HBARC,
                       mom_f/HBARC, phi, self.beam_energy, self.thetalept,
                       self.iform)
        sig = np.array(sig)
        wgt = (
            (var.mom**2 * pke[0] * self.n_z * sig * np.sqrt(MQE**2 + mom_f**2)
             * 2 * np.pi)
            / (var.mom * var.qval)
        )
        return wgt, var

    # TODO: How to best implement this numerically?
    # Currently, this is very unstable from run to run
    def _delta(self, omegap, mom, qval, cost):
        arg = omegap**2 - mom**2 - qval**2 - MQE**2 - 2*mom*qval*cost
        return 1.0/(self.width * np.sqrt(np.pi))*np.exp(-(arg/self.width)**2)

    def _map_vars(self, point):
        """
        Maps a point from the vegas integrator, which takes values in the unit
        ball into the space of physical variables.
        Args:
            point, a point from the integrator

        Returns:
            (variables, ps_wg), where
                * variables is a dict with keys 'mom', 'energy', 'omega', and
                  possibly 'omegap';
                * ps_wgt is the phase space weight; and
        """
        domega = self.wmax - self.wmin
        omega = domega*point[0] + self.wmin
        denergy = self.emax - self.emin
        e_int = denergy*point[1] + self.emin
        dmom = self.pmax - self.pmin
        p_int = dmom*point[2] + self.pmin
        ps_wgt = domega*denergy*dmom
        qval = np.sqrt((2.0 * self.beam_energy * (self.beam_energy - omega)
                        * (1.0 - self.coste)) + omega**2)

        # if FLAGS.folding:
        #     domegap = self.wmax - self.wmin
        #     omegap = domegap*point[3] + self.wmin
        #     ps_wgt = domegap*denergy*dmom
        #     qval_f = np.sqrt(
        #         2.0 * self.beam_energy
        #         * (self.beam_energy - omegap)
        #         * (1.0 - self.coste) + omegap**2)
        #     var = Variables(p_int, e_int, omega, omegap=omegap,
        #                     cost_te=None, qval=qval, qval_f=qval_f)
        #     return var, ps_wgt, domega

        var = Variables(p_int, e_int, omega, omegap=None,
                        cost_te=None, qval=qval, qval_f=None)
        # TODO: Why don't we return domega as the final element of the tuple?
        return var, ps_wgt, None

    def generate_weight(self, point):
        """
        Generates the weight of an event given a point from a vegas integrator.
        The code below intends to follow the following conventions:
        wgt    <--> "weight"             <--> dsigma(omega, E, p)
        ps_wgt <--> "phase space weight" <--> "weight x unit conversion factor"
        wgt_f  <--> "folding weight"     <--> dsigma(omegap, E, p)
        Note that wgt_f should be a function of omegap ("omega prime").
        """
        def _swap_folding(var):
            """Swaps omega <--> omegap and qval <--> qval_f."""
            return Variables(
                var.mom, var.energy, var.omegap, var.omega, var.cost_te,
                var.qval_f, var.qval)

        var, ps_wgt, domega = self._map_vars(point)

        spectral_function = self.pke(var.mom, var.energy)
        wgt, var = self._eval(var, spectral_function)
        wgt *= ps_wgt * TO_NB  # mb -> nb

        # if FLAGS.folding:
        #     swapped_var = _swap_folding(var)
        #     wgt_f, _ = self._eval(swapped_var, spectral_function)
        #     wgt_f *= ps_wgt * TO_NB  # mb -> nb
        #     fold = self.folding(
        #         settings().folding_func,
        #         var.omega,
        #         var.omegap
        #     )
        #     # Note: the term of the form (wgt * transparency) does NOT need a
        #     # factor of domega. In the analytic expression, this term has a
        #     # delta function ~ "delta(omega-omegap) domega". We do this trival
        #     # integral analytically, which eliminates the factor of domega.
        #     wgt_f = (wgt_f * domega * fold) + (wgt * self.folding.transparency)
        #     return wgt_f, var

        return wgt, var

    def generate_momentum(self, var):
        e_out = np.sqrt(var.mom**2 + MQE**2)
        if self.f_g == 1:
            omegat = var.omega
        else:
            omegat = var.omega - var.energy + MQE - e_out
        u_pq = 0.0
        # if FLAGS.folding:
        #     tkin_pf = np.sqrt(var.qval**2+MQE**2)-MQE
        #     if self.folding.kin_max < tkin_pf < self.folding.kin_min:
        #         u_pq = self.folding.kinematic(tkin_pf)

        cost_te = (((omegat+e_out-u_pq)**2 - var.mom**2 - var.qval**2
                    - MQE**2)
                   / (2*var.mom*var.qval))
        if abs(cost_te) > 1:
            logging.debug('omegat = %e, ep = %e, p = %e, qval = %e, mqe = %e'
                          % (omegat, e_out, var.mom, var.qval, MQE))
            return None
        var = _update_variables(var, 'cost_te', cost_te)
        mom_f = np.sqrt(var.mom**2 + var.qval**2
                        + 2 * var.qval * var.mom * var.cost_te)
        epf = np.sqrt(MQE**2 + mom_f**2)
        phi = 2.0*np.pi*np.random.random()

        x_q = var.qval/HBARC
        x_k = var.mom/HBARC
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
