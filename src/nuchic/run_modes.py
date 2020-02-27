""" Module containing all possible run modes.

This module implements a self-registering factory based on the answers found at:
https://stackoverflow.com/questions/55973284/how-to-create-self-registering-factory-in-python


"""

import numpy as np
import seaborn as sns
import pylab as plt

import vectors
import particle
import nucleus
import interactions
import cascade
import logger

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
        density_name = settings().get_run('nuclear_density')
        logger.info(f"Using density '{density_name}'.")
        logger.info(f"Found settings {settings().nucleus_config}.")
        density = densities.NuclearDensity(density_name,
                                           **settings().nucleus_config)

        name = settings().get_run('nucleus')
        binding_energy = settings().nucleus_binding(name)
        fermi_momentum = settings().nucleus_kf(name)
        logger.info(
            f"Building nucleus '{name}' "
            f"with Fermi momentum {fermi_momentum} MeV and "
            f"binding_energy {binding_energy} MeV.")
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
        test_part = particle.Particle(self.pid, momentum, position, -2)
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
        logger.info('Welcome to mean free path')
        # Remove the name argument
        args = args[1:]

        # Delete unused kwargs
        del kwargs

        # Initialize base class and additional variables
        super().__init__()
        self.radius = 10
        self.pid = 2212
        name = settings().get_param('interaction')
        fname = settings().get_param('interaction_file')
        logger.info(f"Creating interaction: '{name}' "
                    f"using interaction data from file: {fname}.")
        interaction = interactions.Interactions.create(name, fname)
        self.fsi = cascade.Cascade(interaction)
        logger.info("Nucleus contains "
                    f"A={self.nucleus.n_nucleons()} total nucleons, "
                    f"Z={self.nucleus.n_protons()} protons, and "
                    f"(A-Z)={self.nucleus.n_neutrons()} neutrons")
        logger.info(f"Binding energy E={self.nucleus.binding_energy()} MeV")
        logger.info(f"Fermi momentum kf={self.nucleus.fermi_momentum()} MeV")
        logger.info(f"Potential energy V={self.nucleus.potential_energy():.2f} MeV")
        logger.info(f"Radius r={self.nucleus.radius():.2f} fm")


    def generate_one_event(self):
        """ Generate one pN or nC event. """

        # Kick in a random direction
        # x in [0,1] --> cos(theta) in [-1,1] --> theta in [0,pi]
        # x in [0,1] --> phi in [0,2*pi]
        x = np.random.random(2)
        theta = np.arccos(2 * x[0] - 1)
        phi = 2 * np.pi * x[1]
        p_kick = settings().beam_energy
        particles = self.nucleus.generate_config()

        # Select a random particle to kick
        kicked_idx = np.random.choice(np.arange(len(particles)))
        self.fsi.set_kicked(kicked_idx)
        kicked_particle = particles[kicked_idx]
        kicked_particle.set_status(-3)
        kicked_particle.set_momentum(vectors.Vector4(
            p_kick * np.sin(theta) * np.cos(phi),
            p_kick * np.sin(theta) * np.sin(phi),
            p_kick * np.cos(theta),
            np.sqrt(kicked_particle.mass()**2.0 + p_kick**2.0)))

        particles = self.fsi.mean_free_path(
            particles,
            self.nucleus.fermi_momentum(),
            self.nucleus.radius()**2)

        return particles

    def finalize(self, events):
        """ Plot a histogram of distance traveled """
        distance_traveled = []
        nhits = 0
        for event in events:
            for aparticle in event:
                if aparticle.status() == -3:
                    distance_traveled.append(aparticle.get_distance_traveled())
                    nhits = nhits + 1
        logger.info(f"nhits / nevents : {nhits} / {len(events)}")
        _, ax = plt.subplots(1)
        sns.distplot(distance_traveled, ax=ax, kde=False)
        ax.set_xlabel(r'Distance traveled [fm]')
        ax.set_ylabel("Counts")
        ax.set_title(
            f"Mean Free Path\n nhits / nevents : {nhits} / {len(events)}"
        )
        plt.show()


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
