""" Implement tests for the Particle class. """

import pytest
import numpy as np

from vectors import Vector4, Vector3
from particle import Particle
# from nuchic.constants import MQE as mN

mN = 938


def test_particle_init():
    """ Test particle initialization. """
    part = Particle(2212, Vector4(4, 3, 2, 1), Vector3(0, 0, 0), 0, [], [])
    assert part.pid() == 2212
    assert part.momentum() == Vector4(4, 3, 2, 1)
    assert part.position() == Vector3(0, 0, 0)
    assert part.status() == 0
    assert part.mothers() == []
    assert part.daughters() == []
    with pytest.raises(Exception):
        Particle(2212, Vector4(4, 3, 2, 1), Vector3(0, 0, 0), None, None, None)


def test_particle_status():
    """ Test particle status codes. """
    part = Particle(pid=2212, momentum=Vector4(4, 3, 2, 1), status=1)
    assert part.is_final()
    assert not part.is_background()
    assert not part.is_propagating()
    part.set_status(0)
    assert not part.is_final()
    assert part.is_background()
    assert not part.is_propagating()
    part.set_status(-1)
    assert not part.is_final()
    assert not part.is_background()
    assert part.is_propagating()


def test_particle_mom():
    """ Test particle momentums. """
    part = Particle(2212, Vector4(3, 2, 1, 4))
    assert part.energy() == 4
    assert part.px() == 3
    assert part.py() == 2
    assert part.pz() == 1
    assert part.mass() == np.sqrt(2)
    assert part.beta() == Vector3(3.0/4.0, 0.5, 0.25)


def test_particle_prop():
    """ Test propagation of particles. """
    part = Particle(2212, Vector4(4, 3, 2, 1), Vector3(0, 0, 0), 0)
    part.propagate(0.1)
    part.back_propagate(0.1)
    assert part.radius() == 0


def test_particle_fzone():
    """ Test particle formation zone. """
    part = Particle(2212, Vector4(np.sqrt(10**2+mN**2), 0, 0, 10),
                    Vector3(0, 0, 0), 0)
    rap = 1
    p_t = 3
    m_t = np.sqrt(p_t**2 + mN**2)
    part_f = Particle(2212, Vector4(m_t*np.cosh(rap), p_t, 0, m_t*np.sinh(rap)),
                      Vector3(0, 0, 0), 0)
    assert not part.is_in_formation_zone()
    part.set_formation_zone(part.momentum(), part_f.momentum())
    f_zone = part.energy()/np.abs(mN**2 - part.momentum()*part_f.momentum())
    print(mN)
    assert part.formation_zone() == f_zone
    assert part.is_in_formation_zone()
