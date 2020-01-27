""" Unit tests for the Nucleus class. """

import numpy as np
import pytest
from nucleus import Nucleus
from particle import Particle


def density():
    """ Dummy density function. """
    protons = [Particle(2212)]*6
    neutrons = [Particle(2112)]*6
    protons.extend(neutrons)
    return protons


def test_nucleus_init():
    """ Initialize nucleus tests. """
    Nucleus(6, 12, 10, 225, density)

    with pytest.raises(Exception):
        Nucleus(12, 6, 10, 225)


def test_nucleus_radius():
    """ Test nucleus radius calculation. """
    nuc = Nucleus(6, 12, 10, 225)
    assert nuc.radius() > 0


def test_nucleus_config():
    """ Test configurations are valid. """
    Z = 6
    A = 12

    nuc = Nucleus(Z, A, 10, 225, density)
    nucleons = nuc.generate_config()
    protons = nuc.protons()
    neutrons = nuc.neutrons()
    assert len(nucleons) == A
    assert len(protons) == Z and len(neutrons) == A-Z


def test_nucleus_momentum():
    """ Test generated momentum are valid. """
    Z = 6
    A = 12
    nuc = Nucleus(Z, A, 10, 225)
    momentum = nuc.generate_momentum()
    assert len(momentum) == 3
    assert np.dot(momentum, momentum) < nuc.fermi_momentum()**2
