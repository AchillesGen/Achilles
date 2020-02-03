#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Implements the Cascade calculation """

import numpy as np
# from absl import logging

from .particle import Particle
from .data.parse_data import GeantData
from .interaction import sigma_pp, sigma_np

from .constants import MQE as mN, GEV, HBARC, MB
from .utils import to_cartesian, timing
from .four_vector import Vec4
from .three_vector import Vec3
from .config import settings


class FSI:
    """ Notes to self:
    * Initialize nucleus
      - Config of p,n
      - Has as input the “kick” (which particle, how much)
    * “Main cascade” functionality
        ** Reabsorption
        ** Exiting
        Notes
        -----
        TODO: Implement binding energy / nuclear potential
        TODO: Implement realistic:
              - initial state configurations (position and momentum)
              - cross sections (distinguishing n from p)
              - reabsorption routine
    """

    def __init__(self, distance, final=True):
        """
        Generates nucleus configuration and kicked nucleon.

        Args:
            energy_transfer: float, energy transfered to the nucleus.
            final: bool, check if nucleon is inside or outside nucleus.
        """
        self.time_step = None
        self.distance = distance
        self.final = final

        # Generate p,n position distribution
        protons, neutrons = settings().nucleus.generate_config()

        # Define particles from configuration
        # NOTE: their 4-momentum is not physical right now
        dummy_mom = Vec4(mN, np.nan, np.nan, np.nan)
        self.nucleons = [
            Particle(
                pid=2212, mom=dummy_mom, pos=Vec3(*x_j)
            ) for x_j in protons
        ]
        self.nucleons += [Particle(pid=2112, mom=dummy_mom,
                                   pos=Vec3(*x_j)) for x_j in neutrons]

        # Cylinder parameters
        self.cylinder_pt1 = 0
        self.cylinder_pt2 = 0

        self.kicked_idxs = []
        self.scatter = False

        self.interactions = GeantData()

    @property
    def number_nucleons(self):
        """ Returns the number of nucleons in the nucleus """
        return len(self.nucleons)

    @property
    def number_protons(self):
        """ Returns the number of protons in the nucleus """
        count = 0
        for particle in self.nucleons:
            if particle.pid == 2212:
                count += 1
        return count

    @property
    def number_neutrons(self):
        """ Returns the number of neutrons in the nucleus """
        count = 0
        for particle in self.nucleons:
            if particle.pid == 2112:
                count += 1
        return count

    @property
    def proton_indices(self):
        """ Returns the indices for the protons """
        indices = []
        for i, particle in enumerate(self.nucleons):
            if particle.pid == 2212:
                indices.append(i)
        return indices

    @property
    def neutron_indices(self):
        """ Returns the indices for the neutrons """
        indices = []
        for i, particle in enumerate(self.nucleons):
            if particle.pid == 2112:
                indices.append(i)
        return indices

    def kick(self, energy_transfer, dsigma):
        '''
        Randomize kicked particle

        Args:
            energy_transfer: float, energy transfered to the nucleus
            dsigma (list): Cross-section for p and n separately
        '''
        self.kicked_idxs = []

        if np.random.random() < dsigma[0]/np.sum(dsigma):
            indices = self.proton_indices
        else:
            indices = self.neutron_indices

        self.kicked_idxs.append(np.random.choice(indices))
        self.nucleons[self.kicked_idxs[0]].status = -1  # propagating nucleon
        self.nucleons[self.kicked_idxs[0]].mom = energy_transfer

    def reset(self):
        '''
        Reset the FSI parameters to begin the next calculation
        '''
        # Generate p,n position distribution
        protons, neutrons = settings().nucleus.generate_config()

        # Define particles from configuration
        # NOTE: their 4-momentum is not physical right now
        dummy_mom = Vec4(mN, np.nan, np.nan, np.nan)
        self.nucleons = [
            Particle(
                pid=2212, mom=dummy_mom, pos=Vec3(*x_j)
            ) for x_j in protons
        ]
        self.nucleons += [Particle(pid=2112, mom=dummy_mom,
                                   pos=Vec3(*x_j)) for x_j in neutrons]

        self.scatter = False
        # Keep outgoing particles after cascade
        # self.outgoing_particles = []

        # Cylinder parameters
        self.cylinder_pt1 = 0
        self.cylinder_pt2 = 0

    def adaptive_step(self, distance):
        '''Adapt the time step so that highest boost particle travels a
        distance `distance'
        '''
        beta = 0
        for idx in self.kicked_idxs:
            if self.nucleons[idx].beta > beta:
                beta = self.nucleons[idx].beta
        self.time_step = distance/(beta*HBARC)  # This is the adapted time step

    @timing
    def __call__(self, max_steps=10000):
        ''' Performs the full propagation of the kicked nucleons inside
        the nucleus. Updates the list of outgoing_particles with
        all status=+1 particles
        '''
        # Overestimate cross section
        sigma = 1000*MB  # 1 barn xsec = 100 fm^2
        # positions = []
        # positions_temp = []
        for step in range(max_steps):
            logging.debug('*******  STEP {} *******'.format(step))
            if self.kicked_idxs == []:
                logging.debug('No more particles propagating - DONE!')
                break
            # Adapt time step
            self.adaptive_step(self.distance)
            # Update formation zones
            for i in self.kicked_idxs:
                if self.nucleons[i].is_in_formation_zone():
                    self.nucleons[i].formation_zone -= self.time_step
            # copy to avoid changing during iteration
            new_kicked_idxs = list(self.kicked_idxs)
            for kick_idx in self.kicked_idxs:
                logging.debug('kick_idx = {}'.format(kick_idx))
                new_kick_idxs = self.interacted(kick_idx, sigma)

                logging.debug('new_kick_idxs = {}'.format(new_kick_idxs))
                new_kick_idxs = new_kick_idxs.tolist()

                if new_kick_idxs:
                    logging.debug('Hit?')
                    (really_did_hit, self.nucleons[kick_idx],
                     self.nucleons[new_kick_idxs[0]]) = \
                        self.generate_final_phase_space(
                            self.nucleons[kick_idx],
                            self.nucleons[new_kick_idxs[0]]
                        )
                    # if it really hit, add index to new kicked index
                    # list and delete duplicates
                    if really_did_hit:
                        logging.debug('Hit!!!!')
                        new_kicked_idxs.append(new_kick_idxs[0])
                        new_kicked_idxs = list(
                            set(new_kicked_idxs))  # Remove duplicates
                        self.scatter = True
                        if not self.final:
                            return

            self.kicked_idxs = new_kicked_idxs
            logging.debug('kicked_idxs = {}'.format(self.kicked_idxs))

            # After-hit checks
            not_propagating = []
            for i, kick_idx in enumerate(self.kicked_idxs):
                # Nucleon becomes final particle if outside nucleus
                if((self.nucleons[kick_idx].pos.mag > settings().nucleus.radius
                   and self.nucleons[kick_idx].status != -2)
                       or self.nucleons[kick_idx].pos.z > settings().nucleus.radius):
                    not_propagating.append(i)
                    if settings().nucleus.escape(self.nucleons[kick_idx]):
                        self.nucleons[kick_idx].status = 1
                        logging.debug('nucleon {} is OOOOOOUT! '
                                      'status: {}'.format(
                                          kick_idx,
                                          self.nucleons[kick_idx].status)
                                      )
                    else:
                        self.nucleons[kick_idx].status = 2
                        logging.debug('nucleon {} is captured! '
                                      'status: {}'.format(
                                          kick_idx,
                                          self.nucleons[kick_idx].status)
                                      )
            # Delete indices of non-propagating particles.
            # Delete in reverse order to avoid shifting elements.
            for i in sorted(not_propagating, reverse=True):
                del self.kicked_idxs[i]

        logging.debug('Number of steps: {}'.format(step))
        stat_list = [n.status for n in self.nucleons]
        logging.debug('Number of final state nucleons: '
                      '{}'.format(sum(stat_list)))
        if -1 in stat_list and self.final:
            for p in self.nucleons:
                print(p.status, p.mom, p.pos)
            logging.fatal(
                "Cascade Failed at step: {}, "
                "has at least one propagating nucleon still".format(step))

        for part in self.nucleons:
            if part.status == -2:
                part.status = 1

        return [n for n in self.nucleons if n.is_final()]

    @staticmethod
    @timing
    def points_in_cylinder(pt1, pt2, radius, position):
        ''' Check if a point is within a cylinder
        Args:
            pt1: initial position vector
            pt2: final position vector
            radius: radius of cylinder
            position: position vector of particle in question
        '''
        pt1 = np.asarray(pt1)
        pt2 = np.asarray(pt2)
        position = np.asarray(position)
        vec = pt2 - pt1
        const = radius * np.linalg.norm(vec)
        return np.logical_and(np.dot(position - pt1, vec) >= 0,
                              np.logical_and(np.dot(position - pt2, vec) <= 0,
                                             np.linalg.norm(
                                                 np.cross(position - pt1, vec),
                                                 axis=-1) <= const))

    @timing
    def interacted(self, idx, sigma):
        ''' Decides if an interaction occurred for a propagating particle
        within given time step

            Parameters
            ----------
                idx : int
                    Index of propagating particle.
                sigma : float
                    Nucleon-nucleon scattering cross section in fm^2

            Returns
            -------
            (key , i), with
            key : bool
                Whether the interaction occured
            i : int
                Particle label

            Notes
            -----
            TODO: This needs to be improved to take into account different
                  cross sections for different nucleai). Maybe this can be
                  an overestimate and the real cross section will be checked in
                  generate_final_phase_space?
        '''

        # Ensure not in formation zone
        if self.nucleons[idx].is_in_formation_zone():
            self.nucleons[idx].propagate(self.time_step)
            return np.array([])

        # Builds up cylinder
        self.cylinder_pt1 = self.nucleons[idx].pos
        logging.debug('Before propagate: {}'.format(self.nucleons[idx]))
        self.nucleons[idx].propagate(self.time_step)
        self.cylinder_pt2 = self.nucleons[idx].pos
        cylinder_r = np.sqrt(sigma/np.pi)
#        key = False
        # Check if any particle (except propagating one) is within cylinder
        # Stops when first is found (not closest one)
        # TODO: optimize with self.n_particles?
        idxs = np.arange(len(self.nucleons))
        idxs = idxs[np.where(idxs != idx)]
        in_cylinder = False

        positions = []
        for i in idxs:
            if self.nucleons[i].is_final() or \
                    self.nucleons[i].is_in_formation_zone():
                idxs = idxs[np.where(idxs != i)]
                continue
            logging.debug('Index: {}, Position: '
                          '{}'.format(i, self.nucleons[i].pos.vec))
            positions.append(self.nucleons[i].pos.vec)
        in_cylinder = self.points_in_cylinder(
            self.cylinder_pt1.array,
            self.cylinder_pt2.array,
            cylinder_r,
            positions
        )
        logging.debug('indices = {}, prop = {}'.format(idxs, idx))
        logging.debug('in_cylinder = {}'.format(in_cylinder))
        return idxs[in_cylinder]

    @timing
    def generate_final_phase_space(self, particle1, particle2):
        ''' Generates phase space (isotropic), checks if particle is inside
        cylinder for proper cross section, checks for Pauli blocking
        (if yes, revert to inital state), set formation zones for
        scattered particles

            Parameters
            ----------
                particle1in : Particle
                    Interating particle 1
                particle2in : Particle
                    Interating particle 2

            Returns
            -------
            (really_did_hit, particle1, particle2), with
            really_did_hit : bool
                False if Pauli blocking occurs (interation did not happen)
            particle1 : Particle
                Outgoing particle 1
            particle2 : Particle
                Outgoing particle 2

            Notes
            -----
            Calls initial phase space generation, gives back fully update
            particles with their corresponding status and formation zones
            (or input particles in case of Pauli blocking)
            TODO: Implement inelastic scattering
            TODO: Discriminate pp, np, nn scatterings (phase space)
            TODO: Implement realistic phase space
        '''

        logging.debug('Before:\nPart1 = {}\nPart2 = {}'.format(particle1,
                                                               particle2))
        # Is particle 2 a background particle? If so, we need to
        # generate it's momentum
        if particle2.is_background():
            # Sort background particle 4-momentum
            p_i = Vec3(*settings().nucleus.generate_momentum())
            energy = np.sqrt(mN**2+p_i.mag2)
            p_mu = Vec4(energy, *p_i.vec)
            particle2.mom = p_mu

        # Start generation of final state phase space
        # Boost back to CoM frame
        total_momentum = particle1.mom+particle2.mom
        boost_vec = total_momentum.boost_vector()

        # Particle 4-momentum in CoM frame
        # See PDG2018 Kinematics Eq. 47.6
        ecm = total_momentum.mass
        pcm = particle1.mom.mom*mN/ecm

        # Fully elastic scattering, protons and neutrons
        # are being treated equally
        # TODO: Inelastic scattering

        # Calculate if particles are inside cylinder with proper momentum
        # dependent cros section.
        # First, we need the momentum of particle 1 in the lab frame.

        # Calculate cross section (flavor dependent)
        if particle1.pid == particle2.pid:
            mode = 'pp'
        else:
            mode = 'np'

        try:
            sigma_p_dependent = self.interactions.cross_section(mode,
                                                                pcm / GEV)
        except ValueError:  # Fall back on hard-coded if not in table
            # logging.warn('Center of mass momentum ({:.3f} MeV) not in table. '
            #              'Falling back on less accurate hard-coded '
            #              'values.'.format(pcm))
            lab_particle1_momentum = \
                particle1.mom.boost_back(particle2.mom.boost_vector()).mom
            if mode == 'pp':
                sigma_p_dependent = sigma_pp(lab_particle1_momentum)
            else:
                sigma_p_dependent = sigma_np(lab_particle1_momentum)

        # Check cylinder:
        # (using global cylinder variables defined in interacted)
        cylinder_r = np.sqrt(sigma_p_dependent/np.pi)
        vec = self.cylinder_pt2.array - self.cylinder_pt1.array
        cylinder_r *= np.linalg.norm(vec)
        particle_r = np.linalg.norm(np.cross(
            particle2.pos.array - self.cylinder_pt1.array, vec))
        if not particle_r <= cylinder_r:
            return False, particle1, particle2

        return self.finalize_momentum(mode, particle1, particle2, boost_vec)

    @timing
    def finalize_momentum(self, mode, particle1, particle2, boost_vec):
        """ Finalize generated momentum. """

        # Boost to center of mass frame
        p1_cm = particle1.mom.boost_back(boost_vec)
        p2_cm = particle2.mom.boost_back(boost_vec)
        ecm = (p1_cm + p2_cm).energy
        pcm = particle1.mom.mom*mN/ecm

        # Generate one outgoing momentum
        @timing
        def make_momentum():
            momentum = np.random.random(3)
            momentum[0] = p1_cm.mom
            momentum[1] = np.radians(
                self.interactions.call(mode, pcm / GEV, momentum[1]))
            momentum[2] = momentum[2]*2*np.pi

            # Three momentum in cartesian coordinates:
            momentum = to_cartesian(momentum)

            return momentum

        momentum = make_momentum()

        # Outgoing 4-momenta
        p1_out = Vec4(*[p1_cm.energy, *momentum])
        p2_out = Vec4(*[p1_cm.energy, *(-momentum)])

        # Get q0 in lab frame for formation zone
        # TODO: what happens in inelastic scatterings?
        # q_lab = p2_out.boost_back(p2_cm.boost_vector()) - \
        #     p2_cm.boost_back(p2_cm.boost_vector())

        # Boost to original frame
        p1_out = p1_out.boost(boost_vec)
        p2_out = p2_out.boost(boost_vec)

        # Check for Pauli blocking and return initial particles if it occurred
        really_did_hit = not(self.pauli_blocking(
            p1_out) or self.pauli_blocking(p2_out))
        if really_did_hit:
            # Assign formation zone
            particle1.set_formation_zone(particle1.mom, p1_out)
            particle2.set_formation_zone(particle2.mom, p2_out)

            # Assign momenta to particles
            particle1.mom = p1_out
            particle2.mom = p2_out

            # Hit background nucleon becomes propagating nucleon
            particle2.status = -1

        logging.debug('After:\nPart1 = {}\nPart2 = {}'.format(particle1,
                                                              particle2))
        return really_did_hit, particle1, particle2

    @timing
    def pauli_blocking(self, four_momentum):
        ''' Checks if Pauli blocking occurred

            Parameters
            ----------
                four_momentum: Fourvector
                    Fourvector of particle in question

            Returns
            -------
            True if Pauli blocking occurred

            Notes
            -----
            Right now, this only checks if magnitude of 3-momentum is
            below Fermi motion
        '''
        # See if Pauli blocking occurs for the proposed interaction
        if four_momentum.mom < settings().nucleus.kf:
            return True
        return False
