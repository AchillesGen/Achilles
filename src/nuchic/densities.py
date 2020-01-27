""" Module containing all possible density models.

This module implements a self-registering factory based on the answers found at:
https://stackoverflow.com/questions/55973284/how-to-create-self-registering-factory-in-python


"""

import pandas as pd
import numpy as np
from scipy.spatial.transform import Rotation

import vectors
import particle

from .config import settings
from .utils import make_path, rand_sphere

class NuclearDensity:
    """ General NuclearDensity class that must be inherited from. """

    class Unknown(Exception):
        """ Error to raise for unknown NuclearDensity. """


    @classmethod
    def get_all_subclasses(cls):
        """ Get all the subclasses recursively of this class. """
        for subclass in cls.__subclasses__():
            yield from subclass.get_all_subclasses()
            yield subclass

    @classmethod
    def _get_name(cls, name):
        return name.lower()

    def __new__(cls, name, *args, **kwargs):
        del args
        del kwargs
        name = cls._get_name(name)
        for subclass in cls.get_all_subclasses():
            if subclass.name == name:
                # Using "object" base class methods avoids recursion here.
                return object.__new__(subclass)
        raise NuclearDensity.Unknown('NuclearDensity "{}" is not known.'.format(name))

    def __init__(self):
        """ General nuclear density initialization. """

    def __call__(self):
        """ Generate a nuclear density configuration. """
        raise NotImplementedError


class NuclearConfiguration(NuclearDensity):
    """ Class that generates configurations from a list of configurations. """
    name = 'configuration'

    def __init__(self, name, **kwargs):
        del name
        super().__init__()

        filename = '{}_configs.out.gz'.format(kwargs['config_type'].upper())
        self.density = pd.read_csv(make_path(filename, 'configurations'),
                                   sep=r'\s+',
                                   names=['pid', 'x', 'y', 'z'],
                                   compression='gzip')

    def __call__(self):
        """ Generate a nuclear configuration based on a configuration file. """
        config_index = np.random.randint(
            0, high=len(self.density.index) / 12)
        idx0 = config_index * 12

        db_tmp = self.density.iloc[idx0:idx0 + 12]
        proton_mask = db_tmp['pid'] == 1
        neutron_mask = db_tmp['pid'] == -1
        protons = np.asarray(db_tmp[proton_mask][['x', 'y', 'z']])
        neutrons = np.asarray(db_tmp[neutron_mask][['x', 'y', 'z']])

        # Rotations using Euler angles in the "x-convention"
        angles = np.random.random(3) * 2 * np.pi
        angles[1] /= 2.
        rotation = Rotation.from_euler('zxz', angles)

        protons = rotation.apply(protons)
        neutrons = rotation.apply(neutrons)

        particles = []
        for proton in protons:
            part = particle.Particle(2212,
                                     vectors.Vector4(0, 0, 0, 0),
                                     vectors.Vector3(proton[0], proton[1], proton[2]))
            particles.append(part)

        for neutron in neutrons:
            part = particle.Particle(2112,
                                     vectors.Vector4(0, 0, 0, 0),
                                     vectors.Vector3(neutron[0], neutron[1], neutron[2]))
            particles.append(part)

        return particles


class NuclearConstant(NuclearDensity):
    """ Class that generates configurations from a constant density. """
    name = 'constant'

    def __init__(self, name, **kwargs):
        del name
        super().__init__()

        self.radius = kwargs['radius']
        self.neutrons = kwargs['n_nucleons']
        self.density = self.neutrons / (4.0/3.0 * np.pi * self.radius**3)

    def __call__(self):
        """ Generate a nuclear configuration with a constant density. """
        nucleons = rand_sphere(self.radius, self.neutrons)
        particles = []
        for nucleon in nucleons:
            part = particle.Particle(2112,
                                     vectors.Vector4(),
                                     vectors.Vector3(nucleon[0], nucleon[1], nucleon[2]))
            particles.append(part)

        return particles
