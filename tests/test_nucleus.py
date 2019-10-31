""" Unit tests for the Nucleus class. """

import numpy as np
import pytest
from nuchic.nucleus import Nucleus


def test_nucleus_init():
    """ Initialize nucleus tests. """
    with pytest.raises(Exception):
        Nucleus(12, 6, 10, 225, 'MF')


def test_nucleus_radius():
    """ Test nucleus radius calculation. """
    nuc = Nucleus(6, 12, 10, 225, 'MF')
    assert nuc.radius > 0


def test_nucleus_config():
    """ Test configurations are valid. """
    Z = 6
    A = 12

    # Mean field
    nuc = Nucleus(Z, A, 10, 225, 'MF')
    protons, neutrons = nuc.generate_config()
    assert len(protons) == Z and len(neutrons) == A-Z
    assert len(np.unique(protons, axis=0)) == Z and len(
        np.unique(neutrons, axis=0)) == A-Z

    # Quantum Monte Carlo
    nuc = Nucleus(Z, A, 10, 225, 'QMC')
    protons, neutrons = nuc.generate_config()
    assert len(protons) == Z and len(neutrons) == A-Z
    assert len(np.unique(protons, axis=0)) == Z and len(
        np.unique(neutrons, axis=0)) == A-Z


def test_nucleus_momentum():
    """ Test generated momentum are valid. """
    Z = 6
    A = 12
    nuc = Nucleus(Z, A, 10, 225, 'MF')
    momentum = nuc.generate_momentum()
    assert len(momentum) == 3
    assert np.dot(momentum, momentum) < nuc.kf**2
