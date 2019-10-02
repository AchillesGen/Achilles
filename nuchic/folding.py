""" Implement the FSI Folding function """
import os

from scipy import interpolate
from absl import flags
from absl import logging
import numpy as np

FLAGS = flags.FLAGS

DIR, FILE = os.path.split(__file__)


class Folding:
    """ Class implementing the FSI folding function.

    Folding is done via a convolution of the cross-section
    with a folding function.
    """

    def __init__(self, width=9, shift=-15):
        logging.info('Initializing folding function')
        with open(os.path.join(DIR, 'pke', 'folding.in')) as infile:
            transparency, hwfold, nwfold = infile.readline().split()
            transparency = float(transparency)
            hwfold = float(hwfold)
            nwfold = int(nwfold)
            self.wfold = np.empty(nwfold)
            fold = np.empty(nwfold)
            for i in range(nwfold):
                self.wfold[i], fold[i] = infile.readline().split()
                logging.debug('wfold[%d] = %e, fold[%d] = %e'
                              % (i, self.wfold[i], i, fold[i]))

        logging.info('TA = %e' % transparency)
        self.transparency = 1.0 - 2.0*np.sum(fold)*hwfold
        logging.info('Norm folding = %e' % transparency)

        self.fold = interpolate.interp1d(self.wfold, fold)

        logging.info('Initializing optical potential')
        with open(os.path.join(DIR, 'pke', 'realOP_12C_EDAI.dat')) as infile:
            npot = int(infile.readline())
            self.kin = np.empty(npot)
            pot = np.empty(npot)
            for i in range(npot):
                self.kin[i], pot[i] = infile.readline().split()
                logging.debug('kin[%d] = %e, potential[%d] = %e'
                              % (i, self.kin[i], i, pot[i]))

        self.kinematic = interpolate.interp1d(self.kin, pot)
        self.width = width
        self.shift = shift

    def noemi(self, omega, omegap):
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

    __dict__ = {'noemi': noemi, 'breit_wigner': breit_wigner}

    def __call__(self, name, omega, omegap):
        func = getattr(self, name, None)
        if func is not None:
            return func(omega, omegap)
        raise NotImplementedError('Requested folding function {} is not '
                                  'defined. Possible choices are {}.'.format(
                                      name, vars(self)))

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
