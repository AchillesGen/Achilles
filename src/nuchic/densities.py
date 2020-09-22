""" Module containing all possible density models.

This module implements a self-registering factory based on the answers found at:
https://stackoverflow.com/questions/55973284/how-to-create-self-registering-factory-in-python


"""

import gzip

import pandas as pd
import numpy as np
from scipy.spatial.transform import Rotation

from . import physics

# from .config import settings
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
        self.max_weight, self.density = read_density(make_path(filename, 'configurations'))

    def __call__(self):
        """ Generate a nuclear configuration based on a configuration file. """
        while True:
            config_index = np.random.randint(
                0, high=len(self.density.index) / 12)
            idx0 = config_index * 12
            weight = self.density.iloc[idx0]['weight']
            if weight/self.max_weight > np.random.uniform(0.0, 1.0):
                break

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
            part = physics.Particle(2212,
                                    physics.Vector4(0, 0, 0, 0),
                                    physics.Vector3(proton[0], proton[1], proton[2]))
            particles.append(part)

        for neutron in neutrons:
            part = physics.Particle(2112,
                                    physics.Vector4(0, 0, 0, 0),
                                    physics.Vector3(neutron[0], neutron[1], neutron[2]))
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
            part = physics.Particle(2112,
                                    physics.Vector4(),
                                    physics.Vector3(nucleon[0], nucleon[1], nucleon[2]))
            particles.append(part)

        return particles


def read_density(fname, file_format='new'):
    if file_format == 'new':
        return _read_density_new(fname)
    elif file_format == 'old':
        return (1, _read_density_old(fname))
    else:
        raise ValueError("Must choose new ")


def _read_density_old(fname):
    return pd.read_csv(fname,
                       sep=r'\s+',
                       names=['pid', 'x', 'y', 'z'],
                       compression='gzip')


def _read_density_new(fname, n_nucleons=12, n_cols=4):
    """
    Reads quantum Monte Carlo output txt files. Output files contain a
    single header line, followed by the data lines. The header contains
    three numbers ('n_configs', 'max_weight', 'min_weight'). The data lines
    appear in 14-line stanzas:
        * 12 lines of four columns ('pid', 'x', 'y', z')
        * 1 line with a weight
        * 1 blank line
    Args:
        fname: str, full path the the file
        n_nucleons: int, defaults to Carbon-12
        n_cols: int, number of columns in data stanzas
    """
    with gzip.open(fname) as ifile:
        # grab header: (n_configs, max_weight, min_weight)
        n_configs, max_weight, _ = np.loadtxt(ifile, max_rows=1)
        n_configs = int(n_configs)
        # make room
        configs = np.zeros((n_configs * n_nucleons, n_cols))
        weights = np.zeros(n_configs * n_nucleons)
        # read data
        for idx in range(n_configs):
            start = idx * n_nucleons
            stop = start + n_nucleons
            # blank line between stanzas is skipped automatically
            configs[start:stop] = np.loadtxt(ifile, max_rows=n_nucleons)
            weights[start:stop] = np.loadtxt(ifile, max_rows=1)
    # back to Pandas
    configs = pd.DataFrame(configs, columns=['pid', 'x', 'y', 'z'])
    weights = pd.DataFrame(weights, columns=['weight'])
    return (max_weight, pd.concat([configs, weights], axis=1))
