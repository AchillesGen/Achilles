""" Parse data from GEANT for pp and np interactions. """

import h5py
import numpy as np
from scipy import interpolate, optimize
# from absl import logging

from ..utils import timing, make_path


class GeantData:
    """ Class to contain information about pp and np interactions.

    The results are interpolated from data from the GEANT program.

    The data can be found at:
    https://github.com/Geant4/geant4/blob/master/source/processes/hadronic/models/coherent_elastic/include/G4LEnp.hh
    and
    https://github.com/Geant4/geant4/blob/master/source/processes/hadronic/models/coherent_elastic/include/G4LEpp.hh

    Interpretation of p_cm can be found at:
    https://github.com/Geant4/geant4/blob/master/source/processes/hadronic/models/coherent_elastic/src/G4HadronElastic.cc

    >>> import numpy as np
    >>> g = GeantData()
    >>> np.allclose(g('np', 0.8, g._interp['np'](90, 0.8)), 90)
    True

    """
    def __init__(self):
        self._names = np.array(['pp', 'np'])
        self._pcm = {}
        self._sig_tot = {}
        self._sig = {}
        self._theta = np.linspace(0.5, 179.5, 180)
        with h5py.File(make_path('GeantData.hdf5', 'data'), 'r') as data:
            self._pcm['pp'] = data['pp/pcm'][:]
            self._pcm['np'] = data['np/pcm'][:]
            self._sig_tot['pp'] = data['pp/sigtot'][:]
            self._sig_tot['np'] = data['np/sigtot'][:]
            self._sig['pp'] = data['pp/sig'][:]
            self._sig['np'] = data['np/sig'][:]

        interp_pp = interpolate.interp2d(
            self._theta, self._pcm['pp'], self._sig['pp'])
        interp_np = interpolate.interp2d(
            self._theta, self._pcm['np'], self._sig['np'])
        self._interp = {'pp': interp_pp, 'np': interp_np}
        cross_section_pp = interpolate.interp1d(
            self._pcm['pp'], self._sig_tot['pp'])
        cross_section_np = interpolate.interp1d(
            self._pcm['np'], self._sig_tot['np'])
        self._cross_section = {'pp': cross_section_pp,
                               'np': cross_section_np}

    def __call__(self, mode, energy, rand):
        """ Return the angle that has a probability rand at the given energy.

            Args:
                - mode: Type of interaction (pp or np).
                - energy: Energy of the center of mass of the interaction.
                - rand: Random number used to determine the angle.

            Returns:
                - angle: Angle in degrees of the outgoing nucleons.
        """
        try:
            return optimize.brentq(
                lambda x: self._interp[mode](x, energy) - rand,
                0.5, 179.5, rtol=1e-8)
        except ValueError:
            # Linearly interpolate between 0 and 0.5
            result = 0.5 / self._interp[mode](0.5, energy) * rand
            logging.warn('Random number {:.3e} outside range of '
                         'f(0.5) = {:.3e} and f(179.5) = {:.3e}. '
                         'Returning angle of {:.3e}.'.format(
                             rand,
                             self._interp[mode](0.5, energy)[0],
                             self._interp[mode](179.5, energy)[0],
                             result[0]))
            return result[0]

    @timing
    def call(self, mode, energy, rand):
        """ Return the angle that has a probability rand at the given energy.

            Args:
                - mode: Type of interaction (pp or np).
                - energy: Energy of the center of mass of the interaction.
                - rand: Random number used to determine the angle.

            Returns:
                - angle: Angle in degrees of the outgoing nucleons.
        """
        return self(mode, energy, rand)

    @timing
    def cross_section(self, mode, energy):
        """ Return the total cross-section at a given center of mass energy."""
        return self._cross_section[mode](energy)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    import timeit
#    t = timeit.Timer("g('np',0.8,0.9)",
#                     setup="from __main__ import GeantData; g = GeantData()")
#    print(t.timeit(number=1000)/1000)
#    time = np.array(t.repeat(100, number=1000))/1000.
#    print(np.mean(time), np.std(time, ddof=1))
