import numpy as np
from nuChic.ThreeVector import Vec3

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
    
    def __init__(self, E=0, px=0, py=0, pz=0):
        self.E = E
        self.px = px
        self.py = py
        self.pz = pz

    def __getitem__(self,i):
        if i == 0: return self.E
        if i == 1: return self.px
        if i == 2: return self.py
        if i == 3: return self.pz
        raise Exception('Vec4D')

    def __repr__(self):
        return 'Vec4({0},{1},{2},{3})'.format(self.E, self.px, self.py, self.pz)

    def __str__(self):
        return '({0},{1},{2},{3})'.format(self.E, self.px, self.py, self.pz)
    
    def __add__(self,v):
        if not isinstance(v,Vec4):
            raise Exception('Vec4')
        return Vec4(self.E+v.E, self.px+v.px, self.py+v.py, self.pz+v.pz)

    def __neg__(self):
        return Vec4(-self.E,-self.px,-self.py,-self.pz)

    def __sub__(self,v):
        return self + -v

    def __mul__(self,v):
        if isinstance(v,Vec4):
            return self.dot(v)
        return Vec4(self.E*v,self.px*v,self.py*v,self.pz*v)

    def __rmul__(self,v):
        if isinstance(v,Vec4):
            return self.dot(v)
        return Vec4(self.E*v,self.px*v,self.py*v,self.pz*v)

    def __truediv__(self,v):
        if isinstance(v,Vec4):
            raise Exception('Vec4')
        return (1.0/v)*self

    def __eq__(self,v):
        if isinstance(v,Vec4):
            if(self.E == v.E and
               self.px == v.px and
               self.py == v.py and
               self.pz == v.pz):
                return True
        return False

    def __pow__(self,n):
        if isinstance(n,int):
            if n%2 == 0:
                return (self.dot(self))**(n/2)
        raise Exception('Vec4')

    def dot(self,v):
        """ Computes the dot product with the four-vector v """
        if not isinstance(v,Vec4):
            raise Exception('Vec4')
        return self.E*v.E - (self.px*v.px + self.py*v.py + self.pz*v.pz)

    def M2(self):
        """ Invariant mass squared of the four-vector p^2 = m^2 """
        return self.dot(self)

    def M(self):
        """ Invariant mass of the four-vector """
        return np.sqrt(self.M2())

    def P2(self):
        """ The square of the spatial three-vector """
        return self.px*self.px+self.py*self.py+self.pz*self.pz

    def P(self):
        """ The magnitude of the spatial three-vector """
        return np.sqrt(self.P2())

    def Pt2(self):
        """ The square of component of the 3-vector transverse to pz """
        return self.px*self.px+self.py*self.py

    def Pt(self):
        """ The magnitude of component of the 3-vector transverse to pz """
        return np.sqrt(self.Pt2())

    def Theta(self):
        """ The polar angle theta of spherical coordinates """
        return np.arccos(self.pz/self.P())

    def Phi(self):
        """ The azimuthal angle phi of spherical coordinates """
        phi = np.arctan2(self.py,self.px)
        if phi < 0:
            phi += 2*np.pi
        return phi

    def Cross(self,v):
        """ The spatial cross product (p x v)_i = eps_{ijk} p_j v_k """
        if not isinstance(v,Vec4):
            raise Exception('Vec4')
        return Vec4(0.0,
                    self.py*v.pz - self.pz*v.py,
                    self.pz*v.px - self.px*v.pz,
                    self.px*v.py - self.py*v.px)

    def BoostVector(self):
        return Vec3(self.px/self.E, self.py/self.E, self.pz/self.E)

    def Boost(self,beta):
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
        if not isinstance(beta,Vec3):
            raise Exception('Vec3')

        beta2 = sum(x*x for x in beta.Vec())
        gamma = 1.0/np.sqrt(1.0-beta2)
        betap = beta[0]*self.px+beta[1]*self.py+beta[2]*self.pz
        gamma2 = (gamma-1.0)/beta2 if beta2 > 0 else 0.0

        return Vec4(gamma*(self.E+betap),
            self.px+gamma2*betap*beta[0]+gamma*beta[0]*self.E,
            self.py+gamma2*betap*beta[1]+gamma*beta[1]*self.E,
            self.pz+gamma2*betap*beta[2]+gamma*beta[2]*self.E)

    def BoostBack(self,beta):
        """Boosts the four-vector by minus beta. See the description of Boost
        for more details.
        Args:
            beta: Vec3, the three-vector defining the boost
        Returns:
            Vec4, the boosted four-vector
        """
        beta = - beta
        return self.Boost(beta)
