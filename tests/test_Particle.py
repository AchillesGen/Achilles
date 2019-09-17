from nuChic.particle import Particle
from nuChic.four_vector import Vec4
from nuChic.three_vector import Vec3
import numpy as np
import pytest

def test_Particle_init():
    part = Particle(2212,Vec4(4,3,2,1),Vec3(0,0,0),0,None,None)
    assert part.pid == 2212
    assert part.mom == Vec4(4,3,2,1)
    assert part.pos == Vec3(0,0,0)
    assert part.status == 0
    assert part.mothers == []
    assert part.daughters == []
    with pytest.raises(Exception):
        Particle(2212,Vec4(4,3,2,1),Vec3(0,0,0),None,None,None)

def test_Particle_status():
    part = Particle(2212,Vec4(4,3,2,1),status=1)
    assert part.is_final()
    assert not part.is_background()
    assert not part.is_propagating()
    part.status = 0
    assert not part.is_final()
    assert part.is_background()
    assert not part.is_propagating()
    part.status = -1
    assert not part.is_final()
    assert not part.is_background()
    assert part.is_propagating()

def test_particle_mom():
    part = Particle(2212,Vec4(4,3,2,1),1,0,None)
    assert part.energy == 4
    assert part.p_x == 3
    assert part.p_y == 2
    assert part.p_z == 1
    assert part.mass == np.sqrt(2)
    assert part.vec() == Vec3(3.0/4.0,0.5,0.25)

def test_Particle_prop():
    part = Particle(2212, Vec4(4,3,2,1), Vec3(0,0,0),0)
    part.propagate(0.1)
    part.back_propagate(0.1)
    assert part.radius == 0

def test_Particle_fzone():
    part = Particle(2212, Vec4(4,3,2,1), Vec3(0,0,0),0)
    assert not part.is_in_formation_zone()
    part.set_formation_zone(1.0,10,10)
    assert part.formation_zone == 1.0/(-10+10**2)
    assert part.is_in_formation_zone()
