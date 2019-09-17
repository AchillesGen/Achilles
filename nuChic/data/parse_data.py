""" Parse data from GEANT for pp and np interactions. """

import os

import h5py
import numpy as np
from scipy import interpolate, optimize

DIR, FILE = os.path.split(__file__)


class GeantData:
    """ Class to contain information about pp and np interactions.

    The results are interpolated from data from the GEANT program.

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
        with h5py.File(os.path.join(DIR, 'GeantData.hdf5'), 'r') as data:
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

    def __call__(self, mode, energy, rand):
        """ Return the angle that has a probability rand at the given energy.

            Args:
                - mode: Type of interaction (pp or np).
                - energy: Energy of the center of mass of the interaction.
                - rand: Random number used to determine the angle.

            Returns:
                - angle: Angle in degrees of the outgoing nucleons.
        """
        return optimize.brentq(lambda x: self._interp[mode](x, energy) - rand,
                               0.5, 179.5, rtol=1e-8)

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


if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    import timeit
#    t = timeit.Timer("g('np',0.8,0.9)",
#                     setup="from __main__ import GeantData; g = GeantData()")
#    print(t.timeit(number=1000)/1000)
#    time = np.array(t.repeat(100, number=1000))/1000.
#    print(np.mean(time), np.std(time, ddof=1))
