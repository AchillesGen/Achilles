from nuChic.Cascade import FSI
from nuChic.Nucleus import Nucleus
import numpy as np
import pytest

def test_Cascade_init():
    pass

def test_Nucleus_radius():
    nuc = Nucleus(6,12,10,225)
    assert nuc.size() > 0

def test_Nucleus_config():
    Z = 6
    A = 12
    nuc = Nucleus(Z,A,10,225)
    protons, neutrons = nuc.generate_config()
    assert len(protons) == Z and len(neutrons) == A-Z
    assert np.all(protons[:,0]**2+protons[:,1]**2+protons[:,2]**2 < nuc.size()**2)
    assert np.all(neutrons[:,0]**2+neutrons[:,1]**2+neutrons[:,2]**2 < nuc.size()**2)
    assert len(np.unique(protons,axis=0)) == Z and len(np.unique(neutrons,axis=0)) == A-Z

def test_Nucleus_momentum():
    Z = 6
    A = 12
    nuc = Nucleus(Z,A,10,225)
    momentum = nuc.generate_momentum()
    assert len(momentum) == 3
    assert np.dot(momentum,momentum) < nuc.kf**2
