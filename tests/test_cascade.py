""" Testing cascade code. """

import unittest
from unittest.mock import patch
import numpy as np
# import pytest

from nuchic.cascade import FSI
from nuchic.particle import Particle

from nuchic.constants import MEV, MQE as mN, FM, HBARC
from nuchic.four_vector import Vec4
from nuchic.three_vector import Vec3

from nuchic.config import settings

NPROTONS = 6
PROTONS = np.random.random((NPROTONS, 3))
NNEUTRONS = 6
NEUTRONS = np.random.random((NNEUTRONS, 3))


def _get_momentums(particles, mode):
    momentum = np.empty_like(particles, dtype=np.float32)
    for i, particle in enumerate(particles):
        momentum[i] = particle.mom[mode]

    return momentum


class TestCascade(unittest.TestCase):
    """ Test the Cascade class. """
    @patch('nuchic.nucleus.Nucleus')
    def setUp(self, mock_nucleus):  # pylint: disable=arguments-differ
        """ Setup needed information. """
        self.k_f = 225 * MEV
        self.mock_nucleus = mock_nucleus
        self.mock_nucleus.generate_config.return_value = (PROTONS, NEUTRONS)
        self.mock_nucleus.kf = self.k_f
        settings().nucleus = self.mock_nucleus

    def test_init(self):
        """ Test cascade initialization. """
        cascade = FSI(1.0)
        self.assertEqual(cascade.number_nucleons, NPROTONS + NNEUTRONS)
        self.assertEqual(cascade.number_protons, NPROTONS)
        self.assertEqual(cascade.number_neutrons, NNEUTRONS)
        self.assertTrue(self.mock_nucleus.generate_config.called_once())

    def test_kick(self):
        """ Test cascade kick. """
        cascade = FSI(1.0)

        e_p = mN + 500*MEV
        p_p = np.sqrt(e_p**2-mN**2)
        energy_transfer = Vec4(e_p, 0, 0, p_p)
        cascade.kick(energy_transfer)

        self.assertEqual(len(cascade.kicked_idxs), 1)
        self.assertEqual(cascade.nucleons[cascade.kicked_idxs[0]].status, -1)
        self.assertEqual(cascade.nucleons[cascade.kicked_idxs[0]].mom.energy,
                         e_p)
        self.assertEqual(cascade.nucleons[cascade.kicked_idxs[0]].mom.p_x, 0)
        self.assertEqual(cascade.nucleons[cascade.kicked_idxs[0]].mom.p_y, 0)
        self.assertEqual(cascade.nucleons[cascade.kicked_idxs[0]].mom.p_z, p_p)

    def test_adaptive_step(self):
        """ Test adaptive_step. """
        cascade = FSI(1.0)

        cascade.nucleons[0].pos = Vec3(0, 0, 0)
        cascade.nucleons[0].mom = Vec4(1000, 0, 0, 500)
        cascade.nucleons[0].status = -1
        cascade.kicked_idxs = [0]
        cascade.adaptive_step(1.0)

        self.assertEqual(cascade.time_step, 1.0/(0.5*HBARC))

        cascade.nucleons[0].propagate(cascade.time_step)

        self.assertEqual(cascade.nucleons[0].pos, Vec3(0, 0, 1.0))

    def test_reset(self):
        """ Test cascade reset. """
        cascade = FSI(1.0)

        e_p = mN + 500*MEV
        p_p = np.sqrt(e_p**2-mN**2)
        energy_transfer = Vec4(e_p, 0, 0, p_p)
        cascade.kick(energy_transfer)
        cascade.reset()

        self.assertEqual(cascade.number_nucleons, NPROTONS + NNEUTRONS)
        self.assertEqual(cascade.number_protons, NPROTONS)
        self.assertEqual(cascade.number_neutrons, NNEUTRONS)
        self.assertEqual(cascade.cylinder_pt1, 0)
        self.assertEqual(cascade.cylinder_pt2, 0)
        energies = _get_momentums(cascade.nucleons, 0)
        output = mN*np.ones_like(energies, dtype=np.float32)
        self.assertTrue(np.allclose(energies, output))
        self.assertTrue(np.all(np.isnan(_get_momentums(cascade.nucleons, 1))))
        self.assertTrue(np.all(np.isnan(_get_momentums(cascade.nucleons, 2))))
        self.assertTrue(np.all(np.isnan(_get_momentums(cascade.nucleons, 3))))

    def test_points_in_cylinder(self):
        """ Test cylinder check. """
        self.assertTrue(FSI.points_in_cylinder(pt1=[0, 0, 0],
                                               pt2=[0, 0, 1],
                                               radius=1,
                                               position=[0.2, 0.2, 0.2]))
        self.assertTrue(not FSI.points_in_cylinder(pt1=[0, 0, 0],
                                                   pt2=[0, 0, 1],
                                                   radius=1,
                                                   position=[0, 1.1, 0.5]))

    def test_pauli_blocking(self):
        """ Test Pauli blocking for the cascade. """
        cascade = FSI(1.0)

        fv1 = Vec4(0, 0, 0, 0.99*self.k_f)
        fv2 = Vec4(0, 0, 0, 1.01*self.k_f)
        self.assertTrue(cascade.pauli_blocking(fv1))
        self.assertTrue(not cascade.pauli_blocking(fv2))

    def test_interacted(self):
        """ Test cascade interaction model. """
        cascade = FSI(1.0)

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
            cascade.adaptive_step(0.2*FM)
            # propagating nucleon
            cascade.nucleons[cascade.kicked_idxs[0]].status = -1
            cascade.nucleons[10] = Particle(
                pid=2212, mom=Vec4(mN, 0, 0, 0), pos=position)
            print(cascade.nucleons[0].pos, cascade.nucleons[10].pos)
            in_cylinder = cascade.interacted(0, sigma=1)
    #        assert in_cylinder==tests[i]
            print(in_cylinder,
                  tests[i], cascade.nucleons[0].pos, cascade.nucleons[10].pos)
