#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from nuChic.FourVector import Vec4
from nuChic.ThreeVector import Vec3
from nuChic.Particle import Particle
from nuChic.Nucleus import Nucleus
from nuChic.Constants import hbarc, MeV, GeV, fm, mN
from copy import deepcopy
   
import logging
logging.basicConfig(level=logging.DEBUG)

class FSI(object):
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
        TODO: Implement realistic initial state configurations (position and momentum), cross sections (distinguishing n from p), reabsorption routine
    """
    def __init__(self, nucleus, energy_transfer, dt):
        """
        Generates nucleus configuration and kicked nucleon.
        
        Args:
            nucleus: nuChic.Nucleus
            energy_transfer: float, energy transfered to the nucleus.
            dt: float, time step
        """
        self.nucleus = nucleus
        self.energy_transfer = energy_transfer
        self.dt = dt
        # Generate p,n position distribution
        protons,neutrons = self.nucleus.generate_config()
        # Define particles from configuration (their 4-momentum is not physical right now)
        dummy_mom = Vec4(mN, np.nan, np.nan, np.nan)
        self.nucleons = [Particle(pid=2212, mom=dummy_mom, pos=Vec3(*x_j)) for x_j in protons]
        self.nucleons += [Particle(pid=2112, mom=dummy_mom, pos=Vec3(*x_j)) for x_j in neutrons]
        # Randomize kicked particle
        self.kicked_idxs = []
        self.kicked_idxs.append(np.random.randint(low=0, high=len(self.nucleons)))
        self.nucleons[self.kicked_idxs[0]].status = -1 # propagating nucleon
        Ep = mN+energy_transfer
        pp = np.sqrt(Ep**2-mN**2)
        self.nucleons[self.kicked_idxs[0]].mom = Vec4(Ep, 0, 0, pp)        
        # Keep outgoing particles after cascade
        self.outgoing_particles = []
        
        
        
    def __call__(self):
        ''' Performs the full propagation of the kicked nucleons inside the nucleus. Updates the list of outgoing_particles with all status=+1 particles
            '''
#        print('First kick: ',self.kicked_idxs)

        sigma = 10 # 0.1 barn xsec = 10 fm^2
        positions=[]
        positions_temp=[]
        for step in range(500):
            if self.kicked_idxs==[]:
                logging.debug('No more particles propagating - DONE!')
                break
            logging.debug('*******  STEP ',step,' *******')
            # Update formation zones
            for i in range(len(self.nucleons)):
                if self.nucleons[i].is_in_formation_zone():
                    self.nucleons[i].formation_zone -= self.dt
            new_kicked_idxs=list(self.kicked_idxs) # copy to avoid changing during iteration
            for kick_idx in self.kicked_idxs:
                logging.debug('kick_idx = ',kick_idx)
                did_hit, new_kick_idx = self.interacted(kick_idx, sigma) 
                if did_hit :
                    logging.debug('Hit?')
                    really_did_hit, self.nucleons[kick_idx], self.nucleons[new_kick_idx] = self.generate_final_phase_space(self.nucleons[kick_idx], self.nucleons[new_kick_idx])
                    # if it really hit, add index to new kicked index list and delete duplicates
                    if really_did_hit :
                        logging.debug('Hit!!!!')
                        new_kicked_idxs.append(new_kick_idx)
                        new_kicked_idxs = list(set(new_kicked_idxs)) # Remove duplicates

#                    #Is it a background particle? If so, we need to generate it's momentum
#                    if not self.nucleons[new_kick_idx].is_propagating() :
#                        # Sort background particle 4-momentum
#                        p_i = Vec3(*self.nucleus.generate_momentum())
#                        energy = np.sqrt(mN**2+p_i.P2())
#                        p_mu = Vec4(energy, *p_i.Vec())
#                        self.nucleons[new_kick_idx].mom=p_mu
#                        # Hit background nucleon becomes propagating nucleon
#                        self.nucleons[new_kick_idx].status=-1
#                        # Add it to list of kicked particles
#                        new_kicked_idxs.append(new_kick_idx)
#                    # Generating outgoing phase space
#                    p1_out, p2_out = self.generate_final_phase_space(self.nucleons[kick_idx],
#                                                                     self.nucleons[new_kick_idx])
                    # Check for Pauli blocking
#                    if self.pauli_blocking(p1_out) or self.pauli_blocking(p2_out) :
                        # Pauli blocking occurred, revert to old configuration
                        
                else:
                    logging.debug('No hit')

            self.kicked_idxs=new_kicked_idxs
            logging.debug('kicked_idxs = ',self.kicked_idxs)

            # After-hit checks
            not_propagating = []
        #    for i in range(len(kicked_idxs)):
            for i, kick_idx in enumerate(self.kicked_idxs):
        #        kick_idx = kicked_idxs[i]
                # Nucleon becomes final particle if
                # (1) is outside nucleus or
                if (self.nucleons[kick_idx].pos.P() > self.nucleus.radius):
                    not_propagating.append(i)
                    self.nucleons[kick_idx].status=1
                    logging.debug('nucleon ', kick_idx,' is OOOOOOUT! status: ',
                                  self.nucleons[kick_idx].status)       
                # (2) has kinetic energy below some barrier energy
                elif (self.nucleons[kick_idx].E()-mN < 30*MeV):
                    not_propagating.append(i)
                    self.nucleons[kick_idx].status=0
                    logging.debug('nucleon ', kick_idx,' is reabsorbed! status: ',
                                  self.nucleons[kick_idx].status)
            # Delete indices of non-propagating particles. 
            # Delete in reverse order to avoid shifting elements.
            for i in sorted(not_propagating, reverse=True):
                del self.kicked_idxs[i]

            # Save positions for 3D plot
            for jj in range(len(self.nucleons)) :
                positions_temp.append(self.nucleons[jj].pos.Vec())
            positions.append(positions_temp)
            positions_temp=[]
            
#            stat_list = [n.status for n in self.nucleons]
#            print('All status: ', stat_list)


        #     if did_hit:
        #         print('idxs:  ', kicked_idxs, new_kick_idx )
        #         print('out momenta: ', mom1, mom2)
        #         print('sanity check', mom1-nucleons[kick_idx].mom, mom2-nucleons[new_kick_idx].mom)
#        print('Number of steps: ',step)
        stat_list = [n.status for n in self.nucleons]
#        print('All status: ', stat_list)
#        print('Number of final state nucleons: ',sum(stat_list))
        if -1 in stat_list : 
            logging.warning("SHIT!!!!") 

        # Record outgoing particles
        self.outgoing_particles = [n for n in self.nucleons if n.is_final()]
        
    @staticmethod
    def points_in_cylinder(pt1, pt2, r, q):
        ''' pt1: initial position vector
            pt2: final position vector
            r: radius of cylinder
            q: position vector of particle in question
        '''
        pt1 = np.asarray(pt1)
        pt2 = np.asarray(pt2)
        q   = np.asarray(q)
        vec = pt2 - pt1
        const = r * np.linalg.norm(vec)
        return (np.dot(q - pt1, vec) >= 0 and np.dot(q - pt2, vec) <= 0 
                and np.linalg.norm(np.cross(q - pt1, vec)) <= const)
    
    @staticmethod
    def to_cartesian(coords):
        ''' Takes spherical coordinates [r, theta, phi] and transform to cartesian coordinates [x,y,z]
        '''
        #r, theta, phi = coords
        r = coords[0]
        theta = coords[1]
        phi= coords[2]
        x = r*np.sin(theta)*np.sin(phi)
        y = r*np.sin(theta)*np.cos(phi)
        z = r*np.cos(theta)
        return np.transpose(np.array([x, y, z]))
        
        
    def interacted(self, idx, sigma):
        ''' Decides if an interaction occurred for a propagating particle within given time step

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
            TODO: This needs to be improved to take into account different cross sections for different nucleai). Maybe this can be an overestimate and the real cross section will be checked in generate_final_phase_space?
        '''

        # Builds up cylinder
        cylinder_pt1 = self.nucleons[idx].pos
        self.nucleons[idx].propagate(self.dt)
        cylinder_pt2 = self.nucleons[idx].pos
        cylinder_r   = np.sqrt(sigma/np.pi)
#        key = False
        # Check if any particle (except propagating one) is within cylinder
        # Stops when first is found (not closest one)
        idxs = np.arange(len(self.nucleons)) # TODO: optimize with self.n_particles?
        in_cylinder = False
        if self.nucleons[idx].is_in_formation_zone():
            return False, np.nan
        for i in idxs[np.where(idxs!=idx)] :
            if self.nucleons[i].is_final() or self.nucleons[i].is_in_formation_zone():
                continue
            position = self.nucleons[i].pos.Vec()
            in_cylinder = self.points_in_cylinder(cylinder_pt1.array(), cylinder_pt2.array(), cylinder_r, position)
            if in_cylinder : 
                # Found particle in cylinder
                break
        return in_cylinder, i

    def generate_final_phase_space(self, particle1in, particle2in) :
        ''' Generate phase space (isotropic), checks for Pauli blocking (if yes, revert to inital state), set formation zones for scattered particles

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
            Calls initial phase space generation, gives back fully update particles with their corresponding status and formation zones (or input particles in case of Pauli blocking)
            TODO: Implement inelastic scattering
            TODO: Discriminate pp, np, nn scatterings (phase space)
            TODO: Implement realistic phase space
        '''

        # We do not want to change the input particles in case Pauli blocking occurred
        particle1 = deepcopy(particle1in)
        particle2 = deepcopy(particle2in)
        
        # Is particle 2 a background particle? If so, we need to generate it's momentum
        if particle2.is_background() :
            # Sort background particle 4-momentum
            p_i = Vec3(*self.nucleus.generate_momentum())
            energy = np.sqrt(mN**2+p_i.P2())
            p_mu = Vec4(energy, *p_i.Vec())
            particle2.mom=p_mu
            # Hit background nucleon becomes propagating nucleon
            particle2.status=-1

        # Start generation of final state phase space
        # Boost back to CoM frame
        total_momentum = particle1.mom+particle2.mom
        boost_vec = total_momentum.BoostVector()
        cm_momentum = total_momentum.BoostBack(boost_vec)
        # Fully elastic scattering, protons and neutrons are being treated equally
        # TODO: Inelastic scattering

        # The original frame is not the lab frame, since both particles have, in general, non-zero momenta    
        # Particle 4-momentum in CoM frame
        p1_cm = particle1.mom.BoostBack(boost_vec)
        p2_cm = particle2.mom.BoostBack(boost_vec)

        # Fix magnitude, (theta, phi) generated isotropically
        momentum = np.random.random(3)
        momentum[0] = p1_cm.P()
        momentum[1] = np.arccos(2*momentum[1] - 1)
        momentum[2] = momentum[2]*2*np.pi

        # Three momentum in cartesian coordinates:
        momentum = self.to_cartesian(momentum)

        # Outgoing 4-momenta
        p1_out = Vec4(*[p1_cm.E, *momentum])
        p2_out = Vec4(*[p1_cm.E, *(-momentum)])

        # Get q0 in lab frame for formation zone
        # TODO: what happens in inelastic scatterings?
        q_lab = p2_out.BoostBack(p2_cm.BoostVector()) - p2_cm.BoostBack(p2_cm.BoostVector())

        # Boost to lab frame
        p1_out=p1_out.Boost(boost_vec)
        p2_out=p2_out.Boost(boost_vec)

        # Assign momenta to particles
        particle1.mom = p1_out
        particle2.mom = p2_out
        
        # Check for Pauli blocking and return initial particles if it occurred
        really_did_hit = not(self.pauli_blocking(p1_out) or self.pauli_blocking(p2_out))
        if really_did_hit :
            # Assign formation zone
            t = q_lab.M2()
            particle1.set_formation_zone(q_lab.E, t, 0.139)
            particle2.set_formation_zone(q_lab.E, t, 0.139)
            #print("form zone = ",foo, q_lab.E, t)
            return really_did_hit, particle1, particle2
        else :
            return really_did_hit, particle1in, particle2in
#            particle1 = deepcopy(particle1in)
#            particle2 = deepcopy(particle2in)

    

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
            Right now, this only checks if magnitude of 3-momentum is below Fermi motion
        '''
        # See if Pauli blocking occurs for the proposed interaction
        if four_momentum.P() < self.nucleus.kf:
            return True
        return False



