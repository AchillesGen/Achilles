import numpy as np
from nuChic.Particle import Particle
from nuChic.Constants import MeV
import pandas as pd

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

        # FIXME: reduced file size for quicker debugs (100,000 configs)
        # If it is Carbon-12, let's use Noemi's configurations.
        # We read it when defining the nucleus and then we only need to pick a configuration when running the cascade
        # 1 million configurations, no header 
        # index   pid    x    y    z
        if Z==6 and A==12: 
            self.c12Density_db = pd.read_csv("/Users/pmachado/Dropbox/Projects/NuGen/FNALNeuGen/configurations/pos_part_in_v2.out",sep='\s+',names=['index','pid','x','y','z'], nrows=1200000)


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

#    def generate_config(self):
#        def to_cartesian(coords):
#            #r, theta, phi = coords
#            r = coords[:,0]
#            theta = coords[:,1]
#            phi= coords[:,2]
#            x = r*np.sin(theta)*np.sin(phi)
#            y = r*np.sin(theta)*np.cos(phi)
#            z = r*np.cos(theta)
#            return np.transpose(np.array([x, y, z]))
#        
#        protons = np.random.random(self.Z*3)
#        protons = protons.reshape(self.Z,3)
#        protons[:,0] = protons[:,0]*self.radius
#        protons[:,1] = np.arccos(2*protons[:,1] - 1)
#        protons[:,2] = protons[:,2]*2*np.pi
#        
#        neutrons = np.random.random((self.A-self.Z)*3)
#        neutrons = neutrons.reshape((self.A-self.Z),3)
#        neutrons[:,0] = neutrons[:,0]*self.radius
#        neutrons[:,1] = np.arccos(2*neutrons[:,1] - 1)
#        neutrons[:,2] = neutrons[:,2]*2*np.pi
#
#        protons = to_cartesian(protons)
#        neutrons = to_cartesian(neutrons)
#
#        return protons, neutrons

    def generate_config(self):
        """ This reads C-12 configuration files only!!!
        """
        def to_cartesian(coords):
            #r, theta, phi = coords
            r = coords[:,0]
            theta = coords[:,1]
            phi= coords[:,2]
            x = r*np.sin(theta)*np.sin(phi)
            y = r*np.sin(theta)*np.cos(phi)
            z = r*np.cos(theta)
            return np.transpose(np.array([x, y, z]))
        
        if not(self.A==12 and self.Z==6):
            raise Exception('Only C-12 is supported for now!')

        config_index = np.random.randint(1, high=100000)
        i0 = (config_index-1)*12
        
        protons = np.asarray(self.c12Density_db.iloc[i0:i0+6][['x','y','z']])
        neutrons = np.asarray(self.c12Density_db.iloc[i0+6:i0+12][['x','y','z']])

        return protons, neutrons



    def generate_momentum(self):
        def to_cartesian(coords):
            #r, theta, phi = coords
            r = coords[0]
            theta = coords[1]
            phi= coords[2]
            x = r*np.sin(theta)*np.sin(phi)
            y = r*np.sin(theta)*np.cos(phi)
            z = r*np.cos(theta)
            return np.transpose(np.array([x, y, z]))
        
        momentum = np.random.random(3)
        momentum[0] = momentum[0]*self.kf
        momentum[1] = np.arccos(2*momentum[1] - 1)
        momentum[2] = momentum[2]*2*np.pi

        momentum = to_cartesian(momentum)

        return momentum
