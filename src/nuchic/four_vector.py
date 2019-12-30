""" Implement the Four Vector class. """

import numpy as np
from .three_vector import Vec3
from .utils import timing


class Vec4:
    """Relativistic four-vector class.

    The names are written with four-momenta in mind, but the functionality
    works equally well in position space. Functions which compute "transverse"
    quantities or angles use pz as defining the azimuthal direction.

    Attributes
        E: the energy, i.e., the temporal component of the four-vector
        px: x-component of the spatial 3-vector
        py: y-component of the spatial 3-vector
        pz: z-component of the spatial 3-vector
    """

    def __init__(self, e=0, px=0, py=0, pz=0):
        self.energy = e
        self.p_x = px
        self.p_y = py
        self.p_z = pz

    def __getitem__(self, i):
        if i == 0:
            return self.energy
        if i == 1:
            return self.p_x
        if i == 2:
            return self.p_y
        if i == 3:
            return self.p_z
        raise Exception('Vec4D')

    def __repr__(self):
        return 'Vec4({0}, {1}, {2}, {3})'.format(self.energy,
                                                 self.p_x,
                                                 self.p_y,
                                                 self.p_z)

    def __str__(self):
        return '({0}, {1}, {2}, {3})'.format(self.energy,
                                             self.p_x,
                                             self.p_y,
                                             self.p_z)

    def __add__(self, other):
        if not isinstance(other, Vec4):
            raise Exception('Vec4')
        return Vec4(self.energy+other.energy,
                    self.p_x+other.p_x,
                    self.p_y+other.p_y,
                    self.p_z+other.p_z)

    def __neg__(self):
        return Vec4(-self.energy, -self.p_x, -self.p_y, -self.p_z)

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        if isinstance(other, Vec4):
            return self.dot(other)
        return Vec4(self.energy*other, self.p_x*other,
                    self.p_y*other, self.p_z*other)

    def __rmul__(self, other):
        if isinstance(other, Vec4):
            return self.dot(other)
        return Vec4(self.energy*other, self.p_x*other,
                    self.p_y*other, self.p_z*other)

    def __truediv__(self, other):
        if isinstance(other, Vec4):
            raise Exception('Vec4')
        return (1.0/other)*self

    def __eq__(self, other):
        if isinstance(other, Vec4):
            if(self.energy == other.energy and
               self.p_x == other.p_x and
               self.p_y == other.p_y and
               self.p_z == other.p_z):
                return True
        return False

    def __pow__(self, power):
        if isinstance(power, int):
            if power % 2 == 0:
                return (self.dot(self))**(power/2)
        raise Exception('Vec4')

    @property
    def array(self):
        return np.array([self.energy, self.p_x, self.p_y, self.p_z])

    def dot(self, other):
        """ Computes the dot product with the four-vector v """
        if not isinstance(other, Vec4):
            raise Exception('Vec4')
        return (self.energy*other.energy
                - (self.p_x*other.p_x
                   + self.p_y*other.p_y + self.p_z*other.p_z))

    @property
    def mass2(self):
        """ Invariant mass squared of the four-vector p^2 = m^2 """
        return self.dot(self)

    @property
    def mass(self):
        """ Invariant mass of the four-vector """
        return np.sqrt(self.mass2)

    @property
    def mom2(self):
        """ The square of the spatial three-vector """
        return self.p_x*self.p_x+self.p_y*self.p_y+self.p_z*self.p_z

    @property
    def mom(self):
        """ The magnitude of the spatial three-vector """
        return np.sqrt(self.mom2)

    @property
    def p_t2(self):
        """ The square of component of the 3-vector transverse to pz """
        return self.p_x*self.p_x+self.p_y*self.p_y

    @property
    def p_t(self):
        """ The magnitude of component of the 3-vector transverse to pz """
        return np.sqrt(self.p_t2)

    @property
    def theta(self):
        """ The polar angle theta of spherical coordinates """
        return np.arccos(self.p_z/self.mom)

    @property
    def phi(self):
        """ The azimuthal angle phi of spherical coordinates """
        phi = np.arctan2(self.p_y, self.p_x)
        if phi < 0:
            phi += 2*np.pi
        return phi

    @property
    def vec3(self):
        """ Return the three vector momentum. """
        return Vec3(self.p_x, self.p_y, self.p_z)

    def set_vec3(self, vec):
        """ Set the three momentum. """
        if not isinstance(vec, Vec3):
            raise Exception('Vec3')
        self.p_x = vec.x
        self.p_y = vec.y
        self.p_z = vec.z

    def cross(self, other):
        """ The spatial cross product (p x v)_i = eps_{ijk} p_j v_k """
        if not isinstance(other, Vec4):
            raise Exception('Vec4')
        return Vec4(0.0,
                    self.p_y*other.p_z - self.p_z*other.p_y,
                    self.p_z*other.p_x - self.p_x*other.p_z,
                    self.p_x*other.p_y - self.p_y*other.p_x)

    def boost_vector(self):
        """ Return the boost vector that would boost the four vector
        to its rest frame. """
        return Vec3(self.p_x/self.energy,
                    self.p_y/self.energy,
                    self.p_z/self.energy)

    @timing
    def boost(self, beta):
        """Boosts the four-fector along the three-vector beta. A discussion of
        the relevant formulae appears, e.g., around Eq (11.19) in in Sec 11.3
        of J.D. Jackson's "Classical Electrodynamics" (3rd Edition). However,
        this function employs a different sign convention for beta, allegedly
        to agree with the conventions of Root.
        Args:
            beta: Vec3, the three vector defining the boost
        Returns:
            Ve4, the boosted four-vector
        """
        if not isinstance(beta, Vec3):
            raise Exception('Vec3')

        beta2 = sum(x*x for x in beta.vec)
        gamma = 1.0/np.sqrt(1.0-beta2)
        betap = beta[0]*self.p_x+beta[1]*self.p_y+beta[2]*self.p_z
        gamma2 = (gamma-1.0)/beta2 if beta2 > 0 else 0.0

        return Vec4(gamma*(self.energy+betap),
                    self.p_x+gamma2*betap*beta[0]+gamma*beta[0]*self.energy,
                    self.p_y+gamma2*betap*beta[1]+gamma*beta[1]*self.energy,
                    self.p_z+gamma2*betap*beta[2]+gamma*beta[2]*self.energy)

    def boost_back(self, beta):
        """Boosts the four-vector by minus beta. See the description of Boost
        for more details.
        Args:
            beta: Vec3, the three-vector defining the boost
        Returns:
            Vec4, the boosted four-vector
        """
        beta = - beta
        return self.boost(beta)
