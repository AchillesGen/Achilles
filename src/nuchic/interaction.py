""" Calculation of pp and np interactions from NASA.
    References:
        - https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20080014212.pdf
"""

import numpy as np
from .constants import MQE as mN, MEV, MB, GEV

# HZETRN parameters for pp, pn, and nn interactions
HZETRN = {'a': 5.0 * MEV,
          'b': 0.199/np.sqrt(MEV),
          'c': 0.451 * MEV**-0.258,
          'd': 25.0 * MEV,
          'e': 134.0 * MEV,
          'f': 1.187 * MEV**-0.35,
          'g': 0.1 * MEV,
          'h': 0.282 * MEV}

# PDG parameters for pp, pn, and nn interactions
PDG = {'Zpp': 33.45 * MB,
       'Zpn': 35.80 * MB,
       'Y1pp': 42.53 * MB,
       'Y1pn': 40.15 * MB,
       'Y2pp': 33.34 * MB,
       'Y2pn': 30.00 * MB,
       'B': 0.308 * MB,
       's1': 1.0 * GEV**2,
       's0': (5.38 * GEV)**2,
       'n1': 0.458,
       'n2': 0.545}

# JWN parameters for pp, pn, and nn interactions
JWN = {'gamma': 52.5 * MB * GEV**(0.16),
       'alpha': 0.00369 / MEV,
       'beta': 0.00895741 * MEV**(-0.8)}


def sigma_pp(plab):
    """ Calculate the pp and nn cross-section given the momentum in the
    lab frame in MeV. """
    tlab = np.sqrt(plab**2+mN**2)-mN
    result = 0
    if plab < 1.8 * GEV:
        if tlab >= 25 * MEV:
            result = (1+HZETRN['a']/tlab) \
                    * (40+109*np.cos(HZETRN['b']*np.sqrt(tlab))
                       * np.exp(-HZETRN['c']*(tlab-HZETRN['d'])**(0.258)))
        else:
            result = np.exp(6.51*np.exp(-(tlab/HZETRN['e'])**(0.7)))
    elif plab <= 4.7 * GEV:
        result = JWN['gamma']/plab**0.16
    else:
        ecm2 = 2*mN*(mN+np.sqrt(plab**2+mN**2))
        result = PDG['Zpp'] + PDG['B']*np.log(ecm2/PDG['s0'])**2 \
            + PDG['Y1pp']*(PDG['s1']/ecm2)**PDG['n1'] \
            - PDG['Y2pp']*(PDG['s1']/ecm2)**PDG['n2']
    return result


def sigma_np(plab):
    """ Calculate the pn cross-section given the momentum in the
    lab frame in MeV. """
    tlab = np.sqrt(plab**2+mN**2)-mN
    result = 0
    if plab < 0.5 * GEV:
        if tlab >= 0.1 * MEV:
            result = 38 + 12500*np.exp(-HZETRN['f']*(tlab-HZETRN['g'])**0.35)
        else:
            result = 26000 * np.exp(-(tlab/HZETRN['h'])**0.3)
    elif plab <= 2.0 * GEV:
        result = 40 + 10*np.cos(JWN['alpha']*plab - 0.943)\
            * np.exp(-JWN['beta']*plab**0.8+2)
    else:
        ecm2 = 2*mN*(mN+np.sqrt(plab**2+mN**2))
        result = PDG['Zpn'] + PDG['B']*np.log(ecm2/PDG['s0'])**2 \
            + PDG['Y1pn']*(PDG['s1']/ecm2)**PDG['n1'] \
            - PDG['Y2pn']*(PDG['s1']/ecm2)**PDG['n2']
    return result
