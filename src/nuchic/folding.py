""" Implement the FSI Folding function """
from collections import namedtuple

from scipy import interpolate
import numpy as np
# from absl import flags
# from absl import logging

from .utils import make_path

# FLAGS = flags.FLAGS


def read_folding_data():
    """
    Reads folding data from folding.in
    The data in this file has the format
    transparency hwfold nfwold
    wfold_0 fold_0
    wfold_1 fold_1
    ...
    In the first row, transparency is the coefficient TA; hwfold is the number
    of steps used to calculate TA (whatever that means); and, nfold is the
    number of points in the remaining two-column data.
    The column wfold is the difference between energies omega and omega' (omega
    prime). The units of wfold are MeV.
    The colum fold is the value of the folding function. Since the folding
    function is used as a convolution with respect to the variable wfold and
    since it shouldn't change the units, the data in the second column should
    have units of 1/energy = 1/MeV.
    """
    folding_params = namedtuple('FoldingParams',
                                ['transparency', 'hwfold', 'nwfold',
                                 'wfold', 'fold'])

    with open(make_path('folding.in', 'pke')) as infile:
        transparency, hwfold, nwfold = infile.readline().split()
        transparency = float(transparency)
        hwfold = float(hwfold)
        nwfold = int(nwfold)
        wfold = np.empty(nwfold)
        fold = np.empty(nwfold)
        for i in range(nwfold):
            wfold[i], fold[i] = infile.readline().split()
            logging.debug(
                'wfold[%d] = %e, fold[%d] = %e' % (i, wfold[i], i, fold[i]))
        return folding_params(transparency, hwfold, nwfold, wfold, fold)


def read_optical_potential():
    """
    Reads optical potential data from realOP_12C_EDAI.dat.
    The data in this file has the format:
    npot
    kin_0 potential_0
    kin_1 potential_1
    ...
    In the first row, npot specifies the number of lines in the file.
    The following two columns are the data itself. The first colum, kin,
    specifies the kinetic energy (KE) of the incoming nucleon.
    TODO: Check w/ Noemi to see if this is for an incoming or outgoing particle
    The second column is the value of the optical potential.
    Both columns have units of MeV.
    """
    optical_pot = namedtuple('OpticalPotential', ['kinetic_in', 'potential'])
    with open(make_path('realOP_12C_EDAI.dat', 'pke')) as infile:
        npot = int(infile.readline())
        kin = np.empty(npot)
        pot = np.empty(npot)
        for i in range(npot):
            kin[i], pot[i] = infile.readline().split()
            logging.debug(
                'kin[%d] = %e, potential[%d] = %e' % (i, kin[i], i, pot[i]))
    return optical_pot(kin, pot)


class Folding:
    """ Class implementing the FSI folding function.

    Folding is done via a convolution of the cross-section
    with a folding function.
    """

    def __init__(self, width=15, shift=-50):
        logging.info('Initializing folding function')
        folding_params = read_folding_data()
        self.wfold = folding_params.wfold
        self.transparency =\
            1.0 - 2.0*np.sum(folding_params.fold)*folding_params.hwfold

        logging.info('TA = %e' % folding_params.transparency)
        logging.info('Norm folding = %e' % folding_params.transparency)

        self.fold = interpolate.interp1d(self.wfold, folding_params.fold)

        logging.info('Initializing optical potential')
        optical_pot = read_optical_potential()
        self.kin = optical_pot.kinetic_in
        self.kinematic = interpolate.interp1d(self.kin, optical_pot.potential)
        self.width = width
        self.shift = shift

    def particle_spectral(self, omega, omegap):
        """ Preform the folding function.

        Args:
            omega: Final transfered energy
            omegap: Initial transfered energy

        Returns:
            folded result
        """
        if abs(omega-omegap) <= self.wmax and abs(omega-omegap) >= self.wmin:
            return self.fold(abs(omega-omegap))
        return 0

    def breit_wigner(self, omega, omegap):
        """ Preform the folding function with a BW.
        Args:
            omega: Final transfered energy
            omegap: Initial transfered energy

        Returns:
            folded result
        """
        return (self.width/np.pi
                / (self.width**2+(omega-omegap-self.shift)**2)
                * (1-self.transparency))

    folding_funcs = ['particle_spectral', 'breit_wigner']

    def __call__(self, name, omega, omegap):
        func = getattr(self, name.lower(), None)
        if func is not None:
            return func(omega, omegap)
        raise NotImplementedError('Requested folding function {} is not '
                                  'defined. Possible choices are {}.'.format(
                                      name, self.folding_funcs))

    @property
    def wmin(self):
        """ Return minimum energy for folding function. """
        return self.wfold[0]

    @property
    def wmax(self):
        """ Return maximum energy for folding function. """
        return self.wfold[-1]

    @property
    def kin_min(self):
        """ Return minimum kinetic energy for folding potential. """
        return self.kin[0]

    @property
    def kin_max(self):
        """ Return maximum kinetic energy for folding potential. """
        return self.kin[-1]
