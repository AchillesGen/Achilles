import numpy as np
from nuChic.constants import HBARC
from nuChic.four_vector import Vec4
from nuChic.ThreeVector import Vec3

class Particle:
    ''' status : int
            -1, 0, 1 for background, propagating and final state particles
    '''
    def __init__(self,pid=0,mom=Vec4(),pos=Vec3(),status=0,mothers=None,daughters=None):
        self.pid = pid
        self.mom = mom
        self.pos = pos
        if status not in [-1,0,1]:
            raise ValueError("Invalid status, must be -1, 0, or 1.")
        self.status = status
        self.mothers = mothers if mothers is not None else []
        self.daughters = daughters if daughters is not None else []
        self.formation_zone = 0

    def vec(self):
        return Vec3(self.mom[1]/self.mom[0],self.mom[2]/self.mom[0],self.mom[3]/self.mom[0])

    def is_background(self):
        if self.status == 0:
            return True
        return False

    def is_propagating(self):
        if self.status == -1:
            return True
        return False
	
    def is_final(self):
        if self.status == 1:
            return True
        return False

    def Px(self):
        return self.mom.px

    def Py(self):
        return self.mom.py

    def Pz(self):
        return self.mom.pz

    def E(self):
        return self.mom.E

    def M(self):
        return self.mom.M()

    def r(self):
        return self.pos.P()

    def propagate(self,time):
        # Propagation of a particle for a given time step
        # The velocity is in terms of beta, time is in units of GeV^-1
        # To get distance need to multiply by hbar*c
        dist = self.mom.P()/self.mom.E*time*hbarc

        theta = self.vec().Theta()
        phi = self.vec().Phi()

        prop_dist_x = dist*np.sin(theta)*np.cos(phi)
        prop_dist_y = dist*np.sin(theta)*np.sin(phi)
        prop_dist_z = dist*np.cos(theta)

        self.prop_dist = Vec3(prop_dist_x,prop_dist_y,prop_dist_z)
        self.pos += self.prop_dist

        return self.prop_dist

    def back_propagate(self):
        self.pos -= self.prop_dist

    def set_formation_zone(self, q0_lab, t, mu):
        ''' When particles interact there is a coherence region in which the interaction takes place.
        In our treatment, we take the interaction to be a single kick. Therefore, to take into account 
        the ``Formation Zone'' we need to turn off any particle interaction after the first kick for
        some time/distance. The in_formation_zone variable should be used as a countdown variable to
        have an interacting particle again after the formation zone.
        Ref.: L. Stodolski, Formation Zone Description in Multiproduction, August 1975.
        '''
        self.formation_zone = q0_lab/(-t+mu**2) # time in GeV^-1
        return self.formation_zone

    def is_in_formation_zone(self):
        if self.formation_zone > 0 :
            return True
        return False
