import numpy as np
from nuChic.Particle import Particle
from nuChic.Constants import MeV

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
            density_N = lambda r: self.nuclear_density/(self.A-self.Z) if r < self.radius else 0
        self.density_P = density_P
        self.density_N = density_N
        self.potential = (np.sqrt((1000)**2 + self.kf**2) - 1000 + 8) * MeV

    def size(self):
        return self.radius

    def escape(self, particle):
        if particle.pos.P() < self.radius:
            return False
        elif np.sqrt(particle.mom.P()**2+particle.mom.M2()) - particle.mom.M() < self.potential:
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

    def generate_config(self):
        protons = np.random.random(Z*3)
        protons.reshape(Z,3)
        protons[:,0] = protons[:,0]*self.radius
        protons[:,1] = np.arccos(2*protons[:,1] - 1)
        protons[:,2] = protons[:,2]*2*np.pi

        neutrons = np.random.random((A-Z)*3)
        neutrons.reshape((A-Z),3)
        neutrons[:,0] = neutrons[:,0]*self.radius
        neutrons[:,1] = np.arccos(2*neutrons[:,1] - 1)
        neutrons[:,2] = neutrons[:,2]*2*np.pi

        return protons, neutrons
