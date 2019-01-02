import queue
import numpy as np

from nuChic.FourVector import Vec4
from nuChic.ThreeVector import Vec3
from nuChic.Particle import Particle
from nuChic.Nucleus import Nucleus
from nuChic.Interaction import *

Interaction_Dict = \
    {# Pions (Pi^0, Pi^+)
     111: [],
     211: [],
     # rhos (rho^0, rho^+)
     113: [],
     213: [],
     # eta and eta'
     221: [],
     331: [],
     # Kaons (K_L, K_S, K^+)
     130: [],
     310: [],
     321: [],
     # D Mesons (D^+, D^0)
     411: [],
     421: [],
     # proton and neutron
     2212: [[sigma_pp, 2212], [sigma_np, 2112]],
     2112: [[sigma_pp, 2112], [sigma_np, 2212]],
     # Deltas (Delta^++, Delta^+, Delta^0, Delta^-)
     2224: [],
     2214: [],
     2114: [],
     1114: [],
     # Lambda and Sigmas (Lambda, Sigma^+, Sigma^0, Sigma^-)
     3122: [],
     3222: [],
     3212: [],
     3112: [],
     # Cascades (0,-)
     3322: [],
     3312: [],
     # Charmed Baryons (Lambda_c, Sigma_c^++, Sigma_c^+, Sigma_c^0)
     4122: [],
     4222: [],
     4212: [],
     4112: [],
    }

class FSI:
    def __init__(self,nucleus, time_step=1):
        ''' Initialize the FSI class.
    
            Input:
            - nucleus: The definition of the nucleus the interaction occured
                       within. This is defined by the nuclear density.
        '''
        self.nucleus = nucleus
        self.queue = queue.Queue()
        self.time_step = time_step

    def __call__(self,particles):
        ''' This function is the driver for the cascade calculation,
            preforming the selection of the nucleons, and their evolution.
        
            Input: 
            - particles: The nucleons propagating through the nucleus from 
                         initial interaction, defined by their 4-momentums,
                         their positions, and their pids
                         
            Output:
            - escaped: A list of particles that escape the nucleus
        '''
        escaped = []
        self.Z = self.nucleus.Z
        self.N = self.nucleus.A - self.nucleus.Z
        
        # Fill the queue with initial particles
        for particle in particles:
            self.queue.put(particle)
            if particle.pid == 2212:
                self.Z -= 1
            elif particle.pid == 2112:
                self.N -= 1
            
        # Main loop over the queue
        while not self.queue.empty():
            particle = self.queue.get()
            prop_dist = particle.propagate(self.time_step)
            
            # Check if particle is still in nucleus
            if particle.r() > self.nucleus.size():
                # Adjust based on potential
                escape = self.nucleus.escape(particle)
                if escape:
                    escaped.append(particle)
                else:
                    self.nucleus.absorb(particle)
            else:
                particles_interaction = self.interact(particle)
                if particles_interaction is None:
                    self.queue.put(particle)
                    continue
                if self.pauli_blocking(particles_interaction):
                    # If blocked just veto the interaction and not undo the step
                    #particle.back_propagate()
                    self.queue.put(particle)
                    continue
                if particles_interaction is not None:
                    if particle.pid == 2212:
                        self.Z += 1
                    elif particle.pid == 2112:
                        self.N += 1
                    for particle_int in particles_interaction:
                        self.queue.put(particle_int)
                        if particle_int.pid == 2212:
                            self.Z -= 1
                        elif particle_int.pid == 2112:
                            self.N -= 1
        return escaped
                        
    def interact(self,particle):
        # Generate initial momentum for hit particle
        phi = 2*np.pi*np.random.uniform()
        ctheta = 2.0*np.random.uniform()-1.0
        stheta = np.sqrt(1.0-ctheta**2)
        px = self.nucleus.kf * stheta * np.cos(phi) * MeV
        py = self.nucleus.kf * stheta * np.sin(phi) * MeV
        pz = self.nucleus.kf * ctheta * MeV
        E = np.sqrt(mN**2+px**2+py**2+pz**2)
        
        p_in = Vec4(E,px,py,pz)
        
        # Calculate all the possible interaction
        sigma = []
        sigma_rho = []
        for interaction in Interaction_Dict[particle.pid]:
            plab = (particle.mom + p_in).M() * GeV
            sigma.append(interaction[0](plab)*1e-1)
            if interaction[1] == 2212:
                sigma_rho.append(sigma[-1]*self.nucleus.density_P(particle.r())*self.Z/self.nucleus.Z)
            elif interaction[1] == 2112:
                sigma_rho.append(sigma[-1]*self.nucleus.density_P(particle.r())*self.N/(self.nucleus.A-self.nucleus.Z))
            
        sigma_total = sum(sigma)
        sigma_rho_total = sum(sigma_rho)
        prob = sigma_rho_total*particle.mom.P()/particle.mom.E * self.time_step
            
        # Accept-reject on the interaction
        # Then select interaction based on cross-section
        if np.random.uniform() > prob:
            return None
        
        interaction_index = np.random.choice(len(Interaction_Dict[particle.pid]),p=sigma/sigma_total)
        interaction = Interaction_Dict[particle.pid][interaction_index]
        
        particle_in = Particle(interaction[1],Vec4(E,px,py,pz),particle.pos)
        
        # Boost to center of mass frame
        cm = particle.mom + particle_in.mom
        beta = cm.BoostVector()
        p1_cm = particle.mom.Boost(-beta)
        p2_cm = particle_in.mom.Boost(-beta)
        
        # Generate the final state and return
        phi = 2*np.pi*np.random.uniform()
        ctheta = 2.0*np.random.uniform()-1.0
        stheta = np.sqrt(1.0-ctheta**2)
        px =  p1_cm.P() * stheta * np.cos(phi)
        py =  p1_cm.P() * stheta * np.sin(phi)
        pz =  p1_cm.P() * ctheta
        E_1 = np.sqrt(p1_cm.M2()+px**2+py**2+pz**2)
        E_2 = np.sqrt(p2_cm.M2()+px**2+py**2+pz**2)
        
        particle1_out = Particle(particle.pid, Vec4(E_1,px,py,pz).Boost(beta), particle.pos)
        particle2_out = Particle(particle_in.pid, Vec4(E_2,-px,-py,-pz).Boost(beta), particle_in.pos)
        
        return [particle1_out,particle2_out]
        
    def pauli_blocking(self, particles_interaction):
        # See if Pauli blocking occurs for the proposed interaction
        for particle in particles_interaction:
            if particle.mom.P() < self.nucleus.kf * MeV:
                return True
        return False
