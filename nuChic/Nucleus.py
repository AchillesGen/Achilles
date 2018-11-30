import numpy as np
from nuChic.Particle import Particle

class Nucleus:
    def __init__(self,Z,A,binding,kf,density_P=None,density_N=None):
        self.Z = Z
        self.A = A
        self.binding = binding
        self.kf = kf
        self.nuclear_density = 0.16
        self.radius = pow(self.A/(4.0/3.0*np.pi*self.nuclear_density),1.0/3.0)
        if density_P is None:
            density_P = lambda r: self.nuclear_density/self.Z if r < self.radius else 0
        if density_N is None:
            density_N = lambda r: self.nuclear_density/self.Z if r < self.radius else 0
        self.density_P = density_P
        self.density_N = density_N
        self.potential = np.sqrt((self.A*1000)**2 + self.kf**2) - self.A*1000 + 8

    def size(self):
        return self.radius

    def escape(self, particle):
        if particle.pos.P() < self.radius:
            return False
        elif particle.mom.P2()/(2*particle.M()) < self.potential:
            return False
        else:
            theta = particle.mom.Theta()
            phi = particle.mom.Phi()

            particle.mom.px -= self.potential*np.sin(theta)*np.cos(phi)
            particle.mom.py -= self.potential*np.sin(theta)*np.sin(phi)
            particle.mom.pz -= self.potential*np.cos(theta)

            return True

    def absorb(self, particle):
        pass
