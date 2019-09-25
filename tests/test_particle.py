""" Implement tests for the Particle class. """

import pytest
import numpy as np

from nuchic.particle import Particle
from nuchic.four_vector import Vec4
from nuchic.three_vector import Vec3
from nuchic.constants import MQE as mN


def test_particle_init():
    """ Test particle initialization. """
    part = Particle(2212, Vec4(4, 3, 2, 1), Vec3(0, 0, 0), 0, None, None)
    assert part.pid == 2212
    assert part.mom == Vec4(4, 3, 2, 1)
    assert part.pos == Vec3(0, 0, 0)
    assert part.status == 0
    assert part.mothers == []
    assert part.daughters == []
    with pytest.raises(Exception):
        Particle(2212, Vec4(4, 3, 2, 1), Vec3(0, 0, 0), None, None, None)


def test_particle_status():
    """ Test particle status codes. """
    part = Particle(2212, Vec4(4, 3, 2, 1), status=1)
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
    """ Test particle momentums. """
    part = Particle(2212, Vec4(4, 3, 2, 1), 1, 0, None)
    assert part.energy == 4
    assert part.p_x == 3
    assert part.p_y == 2
    assert part.p_z == 1
    assert part.mass == np.sqrt(2)
    assert part.vec() == Vec3(3.0/4.0, 0.5, 0.25)


def test_particle_prop():
    """ Test propagation of particles. """
    part = Particle(2212, Vec4(4, 3, 2, 1), Vec3(0, 0, 0), 0)
    part.propagate(0.1)
    part.back_propagate(0.1)
    assert part.radius == 0


def test_particle_fzone():
    """ Test particle formation zone. """
    part = Particle(2212, Vec4(np.sqrt(10**2+mN**2), 0, 0, 10),
                    Vec3(0, 0, 0), 0)
    rap = 1
    p_t = 3
    m_t = np.sqrt(p_t**2 + mN**2)
    part_f = Particle(2212, Vec4(m_t*np.cosh(rap), p_t, 0, m_t*np.sinh(rap)),
                      Vec3(0, 0, 0), 0)
    assert not part.is_in_formation_zone()
    part.set_formation_zone(part.mom, part_f.mom)
    f_zone = part.energy/np.abs(mN**2 - part.mom*part_f.mom)
    assert part.formation_zone == f_zone
    assert part.is_in_formation_zone()
