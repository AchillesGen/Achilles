from .context import nuChic
from nuChic.Particle import Particle
from nuChic.FourVector import Vec4
import numpy as np
import pytest

def test_Particle_init():
    part = Particle(2212,Vec4(4,3,2,1),1,1,1,None,None)
    assert part.pid == 2212
    assert part.mom == Vec4(4,3,2,1)
    assert part.charge == 1
    assert part.I3 == 1
    assert part.status == 1
    assert part.mothers == None
    assert part.daughters == None

def test_Particle_isFinal():
    part = Particle(2212,Vec4(4,3,2,1),status=1)
    assert part.isFinal()
    part.status = 0
    assert not part.isFinal()

def test_particle_Mom():
    part = Particle(2212,Vec4(4,3,2,1),1,None,None)
    assert part.E() == 4
    assert part.Px() == 3
    assert part.Py() == 2
    assert part.Pz() == 1
    assert part.M() == np.sqrt(2)
