""" Module containing all possible run modes.

This module implements a self-registering factory based on the answers found at:
https://stackoverflow.com/questions/55973284/how-to-create-self-registering-factory-in-python


"""

import numpy as np

import vectors
import particle
import nucleus
import interactions
import cascade

import nuchic.densities as densities
from .config import settings

class RunMode:
    """ General RunMode class that must be inherited from. """

    class Unknown(Exception):
        """ Error to raise for unknown RunMode. """

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
        raise RunMode.Unknown('Run mode "{}" is not known.'.format(name))

    def __init__(self):
        """ General run mode initialization. """
        density = densities.NuclearDensity(settings().get_run('nuclear_density'),
                                           **settings().nucleus_config)

        name = settings().get_run('nucleus')
        binding_energy = settings().nucleus_binding(name)
        fermi_momentum = settings().nucleus_kf(name)
        self.nucleus = nucleus.Nucleus.make_nucleus(name,
                                                    binding_energy,
                                                    fermi_momentum,
                                                    density)

    def generate_one_event(self):
        """ Generate one event according to the run mode. """
        raise NotImplementedError

    def finalize(self, events):
        """ Finalize the events for this run mode. """
        raise NotImplementedError


class CalcCrossSection(RunMode):
    """ Class implementing the calculation of pN and nC cross-sections."""
    name = 'xsec'

    def __init__(self, *args, **kwargs):
        # Remove the name argument
        args = args[1:]

        # Initialize base class and additional variables
        super().__init__()
        self.radius = 10
        self.pid = 2212
        interaction = interactions.Interactions.create(settings().get_param('interaction'),
                                                       settings().get_param('interaction_file'))
        self.fsi = cascade.Cascade(interaction)

    def generate_one_event(self):
        """ Generate one pN or nC event. """
        particles = self.nucleus.generate_config()

        # Generate random position in beam
        while True:
            position = self.radius*(2*np.random.rand(2)-1)
            if np.sum(position**2) < self.radius**2:
                break

        # Add test particle to the rest of them
        position = vectors.Vector3(position[0], position[1], -2.5)
        energy = settings().beam_energy
        nucleon_mass = settings().get_param('mn')
        momentum = vectors.Vector4(0, 0, energy, np.sqrt(energy**2+nucleon_mass**2))
        test_part =  particle.Particle(self.pid, momentum, position, -2)
        particles.append(test_part)
        self.fsi.set_kicked(len(particles)-1)
        particles = self.fsi(particles,
                             self.nucleus.fermi_momentum(),
                             2.5**2)

        return particles

    def finalize(self, events):
        """ Convert the events to a total cross-section. """
        xsec = 0
        for event in events:
            count = 0
            for part in event:
                if part.status() == 1:
                    count += 1
            if count != 0:
                xsec += 1

        uncertainty = np.sqrt(xsec)
        xsec *= np.pi*self.radius**2/settings().nevents*10  # fm^2 to mb
        uncertainty *= np.pi*self.radius**2/settings().nevents*10  # fm^2 to mb
        energy = settings().beam_energy**2/(2000)

        print("E: {}\txsec: {} +/- {} mb".format(energy, xsec, uncertainty))


class CalcMeanFreePath(RunMode):
    """ Class implementing the calculation of pN and nC cross-sections."""
    name = 'mfp'

    def __init__(self, *args, **kwargs):
        # Remove the name argument
        args = args[1:]

        # Delete unused kwargs
        del kwargs

        # Initialize base class and additional variables
        super().__init__()

    def generate_one_event(self):
        """ Generate one pN or nC event. """

    def finalize(self, events):
        """ Convert the events to a total cross-section. """


class CalcTransparency(RunMode):
    """ Class implementing the calculation of pN and nC cross-sections."""
    name = 'transparency'

    def __init__(self, *args, **kwargs):
        # Remove the name argument
        args = args[1:]

        # Delete unused kwargs
        del kwargs

        # Initialize base class and additional variables
        super().__init__()
        self.fsi = cascade.Cascade()

    def generate_one_event(self):
        """ Generate one pN or nC event. """

    def finalize(self, events):
        """ Convert the events to a total cross-section. """


class CalcInteractions(RunMode):
    """ Class implementing the calculation of pN and nC cross-sections."""
    name = 'generate'

    def __init__(self, *args, **kwargs):
        # Remove the name argument
        args = args[1:]

        # Delete unused kwargs
        del kwargs

        # Initialize base class and additional variables
        super().__init__()
        self.fsi = cascade.Cascade()

    def generate_one_event(self):
        """ Generate one pN or nC event. """

    def finalize(self, events):
        """ Convert the events to a total cross-section. """
