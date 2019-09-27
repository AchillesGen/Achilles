""" Implement the Particle class. """

import numpy as np
from .constants import HBARC, MQE as mN
from .four_vector import Vec4
from .three_vector import Vec3
from .utils import timing


class Particle:
    ''' status : int
            -1, 0, 1 for background, propagating and final state particles
    '''

    def __init__(self, pid=0, mom=Vec4(), pos=Vec3(),
                 status=0, mothers=None, daughters=None):
        self.pid = pid
        self.mom = mom
        self.pos = pos
        if status not in [-1, 0, 1]:
            raise ValueError("Invalid status, must be -1, 0, or 1.")
        self.status = status
        self.mothers = mothers if mothers is not None else []
        self.daughters = daughters if daughters is not None else []
        self.formation_zone = 0

    def __str__(self):
        return 'Particle({}, {} , {}, {})'.format(self.pid, self.mom,
                                                  self.pos, self.status)

    def vec(self):
        """ Return the momentum 3 vector. """
        return Vec3(self.mom[1]/self.mom[0],
                    self.mom[2]/self.mom[0],
                    self.mom[3]/self.mom[0])

    def is_background(self):
        """ Check if the particle is a background one. """
        if self.status == 0:
            return True
        return False

    def is_propagating(self):
        """ Check if the particle is a propagating. """
        if self.status == -1:
            return True
        return False

    def is_final(self):
        """ Check if the particle is in the final state. """
        if self.status == 1:
            return True
        return False

    @property
    def p_x(self):
        """ Return the x momentum. """
        return self.mom.p_x

    @property
    def p_y(self):
        """ Return the y momentum. """
        return self.mom.p_y

    @property
    def p_z(self):
        """ Return the z momentum. """
        return self.mom.p_z

    @property
    def energy(self):
        """ Return the energy. """
        return self.mom.energy

    @property
    def mass(self):
        """ Return the particle mass. """
        return self.mom.mass

    @property
    def radius(self):
        """ Return the particle distance from origin. """
        return self.pos.mag

    def propagate(self, time):
        """ Propagation of a particle for a given time step
         The velocity is in terms of beta, time is in units of GeV^-1
         To get distance need to multiply by hbar*c"""
        dist = self.mom.mom/self.mom.energy*time*HBARC

        theta = self.vec().theta
        phi = self.vec().phi

        prop_dist_x = dist*np.sin(theta)*np.cos(phi)
        prop_dist_y = dist*np.sin(theta)*np.sin(phi)
        prop_dist_z = dist*np.cos(theta)

        prop_dist = Vec3(prop_dist_x, prop_dist_y, prop_dist_z)
        self.pos += prop_dist

        return prop_dist

    def back_propagate(self, time):
        """ Revert previous propagation. """
        dist = self.mom.mom/self.mom.energy*time*HBARC

        theta = self.vec().theta
        phi = self.vec().phi

        prop_dist_x = dist*np.sin(theta)*np.cos(phi)
        prop_dist_y = dist*np.sin(theta)*np.sin(phi)
        prop_dist_z = dist*np.cos(theta)

        prop_dist = Vec3(prop_dist_x, prop_dist_y, prop_dist_z)
        self.pos -= prop_dist

    @timing
    def set_formation_zone(self, p_in, p_out):
        ''' When particles interact there is a coherence region in which the
        interaction takes place. In our treatment, we take the interaction to
        be a single kick. Therefore, to take into account the `Formation Zone'
        we need to turn off any particle interaction after the first kick for
        some time/distance. The in_formation_zone variable should be used as
        a countdown variable to have an interacting particle again after the
        formation zone.
        Ref.: L. Stodolsky, Formation Zone Description in Multiproduction, 1975
        Ref.: Phys. Rev. C. 86.015505

        Args:
        ----
            - p_in: Momentum of the particle before the interaction
                    in the lab frame
            - p_out: Momentum of the particle after the interaction
                     in the lab frame
        '''
        # formation zone is a time in units of MeV^-1
        self.formation_zone = p_in.energy/np.abs(mN**2-p_in*p_out)
        return self.formation_zone

    def is_in_formation_zone(self):
        """ Check if particle is in formation zone. """
        if self.formation_zone > 0:
            return True
        return False
