import h5py
import numpy as np
from scipy import interpolate, optimize

import os

DIR, FILE = os.path.split(__file__)

class GeantData:
    def __init__(self):
        self._names = np.array(['pp', 'np'])
        self._pcm = {}
        self._sig_tot = {}
        self._sig = {}
        self._theta = np.linspace(0.5,179.5,180)
        with h5py.File(os.path.join(DIR,'GeantData.hdf5'), 'r') as f:
            self._pcm['pp'] = f['pp/pcm'][:]
            self._pcm['np'] = f['np/pcm'][:]
            self._sig_tot['pp'] = f['pp/sigtot'][:]
            self._sig_tot['np'] = f['np/sigtot'][:]
            self._sig['pp'] = f['pp/sig'][:]
            self._sig['np'] = f['np/sig'][:]

        self._interp_pp = interpolate.interp2d(self._theta, self._pcm['pp'], self._sig['pp'])
        self._interp_np = interpolate.interp2d(self._theta, self._pcm['np'], self._sig['np'])
        self._interp = {'pp': self._interp_pp, 'np': self._interp_np}

    def __call__(self, mode, E, r):
        return optimize.brentq(lambda x : self._interp[mode](x, E) - r, 0.5, 179.5, rtol=1e-8)

if __name__ == '__main__':
    import timeit
    t = timeit.Timer("g('np',0.8,0.9)", setup="from __main__ import GeantData; g = GeantData()")
    print(t.timeit(number=1000)/1000)
    time = np.array(t.repeat(100,number=1000))/1000.
    print(np.mean(time), np.std(time,ddof=1))
