""" Testing cascade code. """

from unittest.mock import patch
import numpy as np
# import pytest

from nuchic.cascade import FSI
from nuchic.constants import MEV, MQE as mN
from nuchic.four_vector import Vec4
from nuchic.three_vector import Vec3
from nuchic.particle import Particle


NPROTONS = 6
PROTONS = np.random.random((NPROTONS, 3))
NNEUTRONS = 6
NEUTRONS = np.random.random((NNEUTRONS, 3))


@patch('nuchic.nucleus.Nucleus')
def test_cascade_init(mock_nucleus):
    """ Test cascade initialization. """
    mock_nucleus.generate_config.return_value = (PROTONS, NEUTRONS)
    cascade = FSI(mock_nucleus, 1.0)
    assert cascade.number_nucleons == NPROTONS + NNEUTRONS
    assert cascade.number_protons == NPROTONS
    assert cascade.number_neutrons == NNEUTRONS
    assert mock_nucleus.generate_config.called_once()


@patch('nuchic.nucleus.Nucleus')
def test_kick(mock_nucleus):
    """ Test cascade kick. """
    mock_nucleus.generate_config.return_value = (PROTONS, NEUTRONS)
    cascade = FSI(mock_nucleus, 1.0)

    e_p = mN + 500*MEV
    p_p = np.sqrt(e_p**2-mN**2)
    energy_transfer = Vec4(e_p, 0, 0, p_p)
    cascade.kick(energy_transfer)

    assert len(cascade.kicked_idxs) == 1
    assert cascade.nucleons[cascade.kicked_idxs[0]].status == -1
    assert cascade.nucleons[cascade.kicked_idxs[0]].mom.energy == e_p
    assert cascade.nucleons[cascade.kicked_idxs[0]].mom.p_x == 0
    assert cascade.nucleons[cascade.kicked_idxs[0]].mom.p_y == 0
    assert cascade.nucleons[cascade.kicked_idxs[0]].mom.p_z == p_p


def _get_momentums(particles, mode):
    momentum = np.empty_like(particles, dtype=np.float32)
    for i, particle in enumerate(particles):
        momentum[i] = particle.mom[mode]

    return momentum


@patch('nuchic.nucleus.Nucleus')
def test_reset(mock_nucleus):
    """ Test cascade reset. """
    mock_nucleus.generate_config.return_value = (PROTONS, NEUTRONS)
    cascade = FSI(mock_nucleus, 1.0)

    e_p = mN + 500*MEV
    p_p = np.sqrt(e_p**2-mN**2)
    energy_transfer = Vec4(e_p, 0, 0, p_p)
    cascade.kick(energy_transfer)
    cascade.reset()

    assert cascade.number_nucleons == NPROTONS + NNEUTRONS
    assert cascade.number_protons == NPROTONS
    assert cascade.number_neutrons == NNEUTRONS
    assert cascade.cylinder_pt1 == 0
    assert cascade.cylinder_pt2 == 0
    energies = _get_momentums(cascade.nucleons, 0)
    output = mN*np.ones_like(energies, dtype=np.float32)
    assert np.allclose(energies, output)
    assert np.all(np.isnan(_get_momentums(cascade.nucleons, 1)))
    assert np.all(np.isnan(_get_momentums(cascade.nucleons, 2)))
    assert np.all(np.isnan(_get_momentums(cascade.nucleons, 3)))


def test_points_in_cylinder():
    """ Test cylinder check. """
    assert FSI.points_in_cylinder(pt1=[0, 0, 0], pt2=[0, 0, 1],
                                  radius=1, position=[0.2, 0.2, 0.2])
    assert not(FSI.points_in_cylinder(pt1=[0, 0, 0], pt2=[0, 0, 1],
                                      radius=1, position=[0, 1.1, 0.5]))


@patch('nuchic.nucleus.Nucleus')
def test_pauli_blocking(mock_nucleus):
    """ Test Pauli blocking for the cascade. """
    k_f = 225 * MEV
    mock_nucleus.generate_config.return_value = (PROTONS, NEUTRONS)
    mock_nucleus.kf = k_f
    cascade = FSI(mock_nucleus, 1.0)

    fv1 = Vec4(0, 0, 0, 0.99*k_f)
    fv2 = Vec4(0, 0, 0, 1.01*k_f)
    assert cascade.pauli_blocking(fv1)
    assert not cascade.pauli_blocking(fv2)


@patch('nuchic.nucleus.Nucleus')
def test_interacted(mock_nucleus):
    """ Test cascade interaction model. """
    mock_nucleus.generate_config.return_value = (PROTONS, NEUTRONS)
    cascade = FSI(mock_nucleus, 1.0)

    e_p = 500 * MEV  # massless proton for easy propagation
    p_p = np.sqrt(e_p**2)

    positions = [Vec3(100, 100, 100), Vec3(0, 0, 1.0/2),
                 Vec3(0, 0, 1.5*1.0), Vec3(0, 0.5, 0.5*1.0),
                 Vec3(0, 1.5, 0.5*1.0)]
    tests = [False, True, False, True, False]
    for i, position in enumerate(positions):
        # Put all other nucleaons everyone far
        cascade.nucleons = [Particle(pid=2212, mom=Vec4(mN, 0, 0, 0),
                                     pos=Vec3(100, 100, 100))
                            for j in range(len(cascade.nucleons))]
        # Kicked particle
        cascade.nucleons[0] = Particle(
            pid=2212, mom=Vec4(e_p, 0, 0, p_p), pos=Vec3(0, 0, 0))
        cascade.kicked_idxs = []
        cascade.kicked_idxs.append(0)
        # propagating nucleon
        cascade.nucleons[cascade.kicked_idxs[0]].status = -1
        cascade.nucleons[10] = Particle(
            pid=2212, mom=Vec4(mN, 0, 0, 0), pos=position)
        print(cascade.nucleons[0].pos, cascade.nucleons[10].pos)
        in_cylinder, _ = cascade.interacted(0, sigma=1)
#        assert in_cylinder==tests[i]
        print(in_cylinder,
              tests[i], cascade.nucleons[0].pos, cascade.nucleons[10].pos)
