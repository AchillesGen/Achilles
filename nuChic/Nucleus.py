import os
import numpy as np
from nuChic.particle import Particle
from nuChic.constants import MEV, MQE as mN
import pandas as pd
from scipy.spatial.transform import Rotation


DIR, FILE = os.path.split(__file__)


def main():
    import timeit
    print(
        timeit.timeit(
            'n.generate_config()',
            setup="""
                from __main__ import Nucleus;
                from nuChic.Constants import MEV;
                n = Nucleus(6,12,50*MeV, 225*MEV)""",
            number=100) /
        100)


class Nucleus:
    """
    Basic Nucleus class.
    TODO: Flesh out docs
    """
    def __init__(self, Z, A, binding, kf, density_P=None, density_N=None):
        if Z > A:
            raise ValueError('Protons <= Total Nucelons')

        self.Z = Z
        self.A = A
        self.binding = binding
        self.kf = kf
        self.nuclear_density = 0.16
        self.radius = pow(self.A / (4.0 / 3.0 * np.pi *
                                    self.nuclear_density), 1.0 / 3.0)
        if density_P is None:
            density_P = self._density_P
        if density_N is None:
            density_N = self._density_N
        self.density_P = density_P
        self.density_N = density_N
        self.potential = (np.sqrt(mN**2 + self.kf**2) - mN + 8*MEV)

        # FIXME: reduced file size for quicker debugs (100,000 configs)
        # If it is Carbon-12, let's use Noemi's configurations.
        # We read it when defining the nucleus and then we only need to pick a
        # configuration when running the cascade
        # 1 million configurations, no header
        # index   pid    x    y    z
        if Z == 6 and A == 12:
            self.c12_density_db = pd.read_csv(
                os.path.join(DIR, 'configurations', 'pos_part_in_v2.out.gz'),
                sep=r'\s+',
                names=['index', 'pid', 'x', 'y', 'z'],
                compression='gzip'
            )
#            self.c12Density_db = pd.read_csv("/Users/pmachado/Dropbox/Projects/NuGen/FNALNeuGen/configurations/pos_part_in_v2.out",sep='\s+',names=['index','pid','x','y','z'], nrows=1200000)

    def _density_P(self, r):
        return self.nuclear_density / self.Z if r < self.radius else 0

    def _density_N(self, r):
        return self.nuclear_density / \
            (self.A - self.Z) if r < self.radius else 0

    def size(self):
        # TODO is this function really necessary?
        # It just seems like an alias for self.radius
        return self.radius

    def escape(self, particle):
        """Check whether or not the particle escaped.
        If the particle escapes, the spatial components of the momentum are
        updated to account for the influence of the potential.
        """
        # Is the position still within the nucleus?
        if particle.pos.p < self.radius:
            return False
        # Is the "kinetic energy" less than the binding potential?
        energy_total = np.sqrt(particle.mom.mom2 + particle.mom.mass2)
        energy_kinetic = energy_total - particle.mom.mass
        if energy_kinetic < self.potential:
            return False
        # Then the particle has escaped. Update the momentum to account for
        # the effect of the potential
        theta = particle.mom.theta
        phi = particle.mom.phi
        particle.mom.p_x -= self.potential * np.sin(theta) * np.cos(phi)
        particle.mom.p_y -= self.potential * np.sin(theta) * np.sin(phi)
        particle.mom.p_z -= self.potential * np.cos(theta)
        return True

    def absorb(self, particle):
        """TODO Implement absorb, add doc"""
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
        """ This reads C-12 configuration files only!!!"""
        def to_cartesian(coords):
            """Convert spherical coordinates to cartesian coordinates.
            TODO This function is not used. Delete?
            """
            # r, theta, phi = coords
            r = coords[:, 0]
            theta = coords[:, 1]
            phi = coords[:, 2]
            x = r * np.sin(theta) * np.sin(phi)
            y = r * np.sin(theta) * np.cos(phi)
            z = r * np.cos(theta)
            return np.transpose(np.array([x, y, z]))

        if not(self.A == 12 and self.Z == 6):
            raise Exception('Only C-12 is supported for now!')

        config_index = np.random.randint(
            0, high=len(self.c12_density_db.index) / 12)
        idx0 = config_index * 12

        db_tmp = self.c12_density_db.iloc[idx0:idx0 + 12]
        proton_mask = db_tmp['pid'] == 1
        neutron_mask = db_tmp['pid'] == 2
        protons = np.asarray(db_tmp[proton_mask][['x', 'y', 'z']])
        neutrons = np.asarray(db_tmp[neutron_mask][['x', 'y', 'z']])

        # Rotations using Euler angles in the "x-convention"
        angles = np.random.random(3) * 2 * np.pi
        angles[1] /= 2.
        rotation = Rotation.from_euler('zxz', angles)

        protons = rotation.apply(protons)
        neutrons = rotation.apply(neutrons)

        return protons, neutrons

    def generate_momentum(self):
        """Generate a random momentum 3-vector in cartesian coordinates."""
        def to_cartesian(coords):
            """Convert spherical coordinates to cartesian coordinates."""
            # r, theta, phi = coords
            r = coords[0]
            theta = coords[1]
            phi = coords[2]
            x = r * np.sin(theta) * np.sin(phi)
            y = r * np.sin(theta) * np.cos(phi)
            z = r * np.cos(theta)
            return np.transpose(np.array([x, y, z]))

        momentum = np.random.random(3)
        # Random numbers fist to spherical coords
        momentum[0] = momentum[0] * self.kf
        momentum[1] = np.arccos(2 * momentum[1] - 1)
        momentum[2] = momentum[2] * 2 * np.pi
        # and then to cartesian
        momentum = to_cartesian(momentum)
        return momentum


if __name__ == '__main__':
    main()
