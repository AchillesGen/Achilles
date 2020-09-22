""" Testing cascade code. """

import unittest
from unittest.mock import patch
import numpy as np
# import pytest

from nuchic.physics import Vector4, Vector3
from nuchic.physics import Particle, ParticleStatus
from nuchic.physics import Nucleus
from nuchic._nuchic import Cascade
from nuchic._nuchic import Interactions

from nuchic._nuchic.utilities.constants import mN, hbarc
from nuchic.config import settings

MEV = 1
FM = 1

NPROTONS = 6
PROTONS = 2.5*np.random.random((NPROTONS, 3))
NNEUTRONS = 6
NEUTRONS = 2.5*np.random.random((NNEUTRONS, 3))


def _get_momentums(particles, mode):
    momentum = np.empty_like(particles, dtype=np.float32)
    for i, particle in enumerate(particles):
        momentum[i] = particle.momentum()[mode]

    return momentum


def density():
    """ Dummy density function. """
    nucleons = []
    for position in PROTONS:
        nucleons.append(Particle(2212, position=Vector3(position)))
    for position in NEUTRONS:
        nucleons.append(Particle(2112, position=Vector3(position)))
    return nucleons


class MockInteractions(Interactions):
    def is_registered(self):
        return True

    def cross_section(self, p1, p2):
        return 10

    def make_momentum(self, mode, x, y, rans):
        return Vector3(10, 10, 10) 


class TestCascade(unittest.TestCase):
    """ Test the Cascade class. """

    def test_init(self):
        """ Test cascade initialization. """
        mock_interact = MockInteractions()
        cascade = Cascade(mock_interact, Cascade.Gaussian)
        self.assertIsInstance(cascade, Cascade)

    def test_kick(self):
        """ Test cascade kick. """
        mock_interact = MockInteractions()
        mock_nucleus = Nucleus(6, 12, 10, 225, "c12.prova.txt",
                               Nucleus.Global, density)
        mock_nucleus.generate_config()
        cascade = Cascade(mock_interact, Cascade.Gaussian)

        e_p = mN + 500*MEV
        p_p = np.sqrt(e_p**2-mN**2)
        energy_transfer = Vector4(0, 0, p_p, e_p)
        particles = mock_nucleus.nucleons()
        cascade.kick(mock_nucleus, energy_transfer, [1, 0])
        kicked_parts = mock_nucleus.nucleons()

        kicked = []
        for i, part in enumerate(kicked_parts):
            if part.status() == ParticleStatus.propagating:
                kicked.append(part)
                idx = i

        final_momentum = energy_transfer + particles[idx].momentum()
        self.assertEqual(len(kicked), 1)
        self.assertEqual(final_momentum, kicked_parts[idx].momentum())

    # This is a private function in C++,
    # maybe add tests in the C++ code to handle these things?
    # def test_adaptive_step(self):
    #     """ Test adaptive_step. """
    #     cascade = Cascade(1.0)

    #     cascade.nucleons[0].pos = Vector3(0, 0, 0)
    #     cascade.nucleons[0].mom = Vector4(1000, 0, 0, 500)
    #     cascade.nucleons[0].status = -1
    #     cascade.kicked_idxs = [0]
    #     cascade.adaptive_step(1.0)

    #     self.assertEqual(cascade.time_step, 1.0/(0.5*HBARC))

    #     cascade.nucleons[0].propagate(cascade.time_step)

    #     self.assertEqual(cascade.nucleons[0].pos, Vector3(0, 0, 1.0))

    # Need to figure out how to test this. It simply ensures the kicked_idxs is size 0
    # def test_reset(self):
    #     """ Test cascade reset. """
    #     cascade = Cascade(1.0)

    #     e_p = mN + 500*MEV
    #     p_p = np.sqrt(e_p**2-mN**2)
    #     energy_transfer = Vector4(e_p, 0, 0, p_p)
    #     cascade.kick(energy_transfer)
    #     cascade.reset()

    #     self.assertEqual(cascade.number_nucleons, NPROTONS + NNEUTRONS)
    #     self.assertEqual(cascade.number_protons, NPROTONS)
    #     self.assertEqual(cascade.number_neutrons, NNEUTRONS)
    #     self.assertEqual(cascade.cylinder_pt1, 0)
    #     self.assertEqual(cascade.cylinder_pt2, 0)
    #     energies = _get_momentums(cascade.nucleons, 0)
    #     output = mN*np.ones_like(energies, dtype=np.float32)
    #     self.assertTrue(np.allclose(energies, output))
    #     self.assertTrue(np.all(np.isnan(_get_momentums(cascade.nucleons, 1))))
    #     self.assertTrue(np.all(np.isnan(_get_momentums(cascade.nucleons, 2))))
    #     self.assertTrue(np.all(np.isnan(_get_momentums(cascade.nucleons, 3))))

    # No longer use cylinder, and plane checking is private and not exposed to python
    # def test_points_in_cylinder(self):
    #     """ Test cylinder check. """
    #     self.assertTrue(Cascade.points_in_cylinder(pt1=[0, 0, 0],
    #                                            pt2=[0, 0, 1],
    #                                            radius=1,
    #                                            position=[0.2, 0.2, 0.2]))
    #     self.assertTrue(not Cascade.points_in_cylinder(pt1=[0, 0, 0],
    #                                                pt2=[0, 0, 1],
    #                                                radius=1,
    #                                                position=[0, 1.1, 0.5]))

    # This is also private in the C++ code
    # def test_pauli_blocking(self):
    #     """ Test Pauli blocking for the cascade. """
    #     cascade = Cascade(1.0)

    #     fv1 = Vector4(0, 0, 0, 0.99*self.k_f)
    #     fv2 = Vector4(0, 0, 0, 1.01*self.k_f)
    #     self.assertTrue(cascade.pauli_blocking(fv1))
    #     self.assertTrue(not cascade.pauli_blocking(fv2))

    # TODO: Add this into its own test suite for ensuring the interactions are correct
    # def test_interacted(self):
    #     """ Test cascade interaction model. """
    #     cascade = Cascade(1.0)

    #     e_p = 500 * MEV  # massless proton for easy propagation
    #     p_p = np.sqrt(e_p**2)

    #     positions = [Vector3(100, 100, 100), Vector3(0, 0, 1.0/2),
    #                  Vector3(0, 0, 1.5*1.0), Vector3(0, 0.5, 0.5*1.0),
    #                  Vector3(0, 1.5, 0.5*1.0)]
    #     tests = [False, True, False, True, False]
    #     for i, position in enumerate(positions):
    #         # Put all other nucleaons everyone far
    #         cascade.nucleons = [Particle(pid=2212, mom=Vector4(mN, 0, 0, 0),
    #                                      pos=Vector3(100, 100, 100))
    #                             for j in range(len(cascade.nucleons))]
    #         # Kicked particle
    #         cascade.nucleons[0] = Particle(
    #             pid=2212, mom=Vector4(e_p, 0, 0, p_p), pos=Vector3(0, 0, 0))
    #         cascade.kicked_idxs = []
    #         cascade.kicked_idxs.append(0)
    #         cascade.adaptive_step(0.2*FM)
    #         # propagating nucleon
    #         cascade.nucleons[cascade.kicked_idxs[0]].status = -1
    #         cascade.nucleons[10] = Particle(
    #             pid=2212, mom=Vector4(mN, 0, 0, 0), pos=position)
    #         print(cascade.nucleons[0].pos, cascade.nucleons[10].pos)
    #         in_cylinder = cascade.interacted(0, sigma=1)
    #         assert in_cylinder==tests[i]
    #         print(in_cylinder,
    #               tests[i], cascade.nucleons[0].pos, cascade.nucleons[10].pos)
