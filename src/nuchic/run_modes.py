""" Module containing all possible run modes.

This module implements a self-registering factory based on the answers found at:
https://stackoverflow.com/questions/55973284/how-to-create-self-registering-factory-in-python


"""

import numpy as np

import seaborn as sns
import pylab as plt
import scipy.stats

import vectors
import particle
import nucleus
import interactions
import cascade
import logger

import nuchic.densities as densities
from .config import settings


class ConstantInteraction(interactions.Interactions):
    """ Interaction model that always returns a constant value. """
    _registered = True

    def __init__(self, xsec):
        interactions.Interactions.__init__(self)
        self.xsec = xsec
        logger.info(f"ConstantInteraction: xsec {self.xsec}")

    def CrossSection(self, part1, part2):  # pylint: disable=unused-argument
        """ Return the cross-section value. """
        return self.xsec

    def IsRegistered(self):  # pylint: disable=invalid-name
        """ Check if the interaction has been registered in the factory. """
        return self._registered

    def MakeMomentum(self, same_pid, p1_cm, pcm, rans):  # pylint: disable=unused-argument
        """ Make the momentum for the outgoing nucleons. """
        return vectors.Vector3()


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
        instance_name = self.__class__.__name__
        density_name = settings().get_run('nuclear_density')
        logger.info(f"{instance_name}: using density '{density_name}'.")
        logger.info(f"{instance_name}: using nuclear configuration settings "
                    f"{settings().nucleus_config}.")
        density = densities.NuclearDensity(density_name,
                                           **settings().nucleus_config)

        name = settings().get_run('nucleus')
        binding_energy = settings().nucleus_binding(name)
        fermi_momentum = settings().nucleus_kf(name)
        density_file = settings().get_param('density_file')
#        try:
        self.fermi_gas= nucleus.Nucleus.__dict__[settings().get_param("fermi_gas")]

#        except KeyError as error:
#            logger.error(f"Invalid FG model {error}")
#            raise
    
        logger.info(
            f"{instance_name}: Building nucleus '{name}' "
            f"with density file {density_file} and "
            f"binding_energy {binding_energy} MeV."
            f"Fermi gas model implemented {self.fermi_gas}")
        self.nucleus = nucleus.Nucleus.make_nucleus(name,
                                                    binding_energy, fermi_momentum,
                                                    density_file, self.fermi_gas,
                                                    density)
        
        try:
            self.cascade_prob = cascade.Cascade.__dict__[settings().get_param("cascade_prob")]
        except KeyError as error:
            logger.error(f"Invalid probability model {error}")
            raise

        logger.info(f"{instance_name}: using cascade probability: '{self.cascade_prob}'.")

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

        # Delete unused kwargs
        del kwargs

        # Initialize base class and additional variables
        super().__init__()
        self.radius = 10
        self.pid = 2212
        interaction = interactions.Interactions.create(settings().get_param('interaction'),
                                                       settings().get_param('interaction_file'))
        self.fsi = cascade.Cascade(interaction, self.cascade_prob)

    def generate_one_event(self):
        """ Generate one pN or nC event. """
        self.nucleus.generate_config() 
        particles = self.nucleus.nucleons()

        # Generate random position in beam
        while True:
            position = self.radius*(2*np.random.rand(2)-1)
            if np.sum(position**2) < self.radius**2:
                break

        # Add test particle to the rest of them
        position = vectors.Vector3(position[0], position[1], -6.5)
        energy = settings().beam_energy
        nucleon_mass = settings().get_param('mn')
        momentum = vectors.Vector4(0, 0, energy, np.sqrt(energy**2+nucleon_mass**2))
        test_part = particle.Particle(self.pid,
                                      momentum,
                                      position,
                                      particle.ParticleStatus.external_test)

        # Set up kicked
        self.fsi.set_kicked(len(particles))
        particles.append(test_part)
        self.nucleus.set_nucleons(particles)
        self.fsi.nuwro(self.nucleus)

        return self.nucleus.nucleons()

    def finalize(self, events):
        """ Convert the events to a total cross-section. """
        xsec = 0
        for event in events:
            count = 0
            for part in event:
                if part.status() == particle.ParticleStatus.escaped:
                    count += 1
            if count != 0:
                xsec += 1

        print(f"nhits: {xsec}")
        uncertainty = np.sqrt(xsec)
        xsec *= np.pi*self.radius**2/settings().nevents*10  # fm^2 to mb
        uncertainty *= np.pi*self.radius**2/settings().nevents*10  # fm^2 to mb
        energy = settings().beam_energy**2/(2000)

        print("E: {}\txsec: {} +/- {} mb".format(energy, xsec, uncertainty))
        return xsec, uncertainty


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
        logger.info("CalcMeanFreePath: overriding nuclear radius "
                    "with value from run card.")
        self.nucleus.set_radius(settings().nucleus_config['radius'])
        self.radius = self.nucleus.radius
        self.pid = 2212
        name = settings().get_param('interaction')
        # fname = str(settings().get('xsec'))
        logger.info(f"CalcMeanFreePath: creating interaction: '{name}'.")
        # interaction = interactions.Interactions.create(name, fname)
        interaction = ConstantInteraction(xsec=settings().get('xsec'))

        self.fsi = cascade.Cascade(interaction, self.cascade_prob)
        del interaction
        logger.info(
            "CalcMeanFreePath: nucleus contains "
            f"A={self.nucleus.n_nucleons()} total nucleons, "
            f"Z={self.nucleus.n_protons()} protons, and "
            f"(A-Z)={self.nucleus.n_neutrons()} neutrons.")
        logger.info(
            "CalcMeanFreePath: nucleus has "
            f"Binding energy E={self.nucleus.binding_energy()} MeV, "
            f"Fermi momentum kf={self.nucleus.fermi_momentum()} MeV, "
            f"Potential energy V={self.nucleus.potential_energy():.2f} MeV, and "
            f"Radius r={self.nucleus.radius():.2f} fm.")

    def generate_one_event(self):
        """ Generate one pN or nC event. """

        # Kick in a random direction
        # x in [0,1] --> cos(theta) in [-1,1] --> theta in [0,pi]
        # x in [0,1] --> phi in [0,2*pi]
        x = np.random.random(2)
        theta = np.arccos(2 * x[0] - 1)
        phi = 2 * np.pi * x[1]
        p_kick = settings().beam_energy
        self.nucleus.generate_config()
        particles = self.nucleus.nucleons()

        # Place a test particle in the center and give it a kick
        position = vectors.Vector3(0.0, 0.0, 0.0)
        nucleon_mass = settings().get_param('mn')
        momentum = vectors.Vector4(
            p_kick * np.sin(theta) * np.cos(phi),
            p_kick * np.sin(theta) * np.sin(phi),
            p_kick * np.cos(theta),
            np.sqrt(p_kick**2.0+nucleon_mass**2))
        test_part = particle.Particle(2212,
                                      momentum,
                                      position,
                                      particle.ParticleStatus.internal_test)

        # Evolve the nucleus
        self.fsi.set_kicked(len(particles))
        particles.append(test_part)
        self.nucleus.set_nucleons(particles)
        self.fsi.mean_free_path(self.nucleus)

        return self.nucleus.nucleons()

    def finalize(self, events):
        """ Plot a histogram of distance traveled """
        distance_traveled = []
        nhits = 0
        for event in events:
            for aparticle in event:
                if aparticle.status() == particle.ParticleStatus.internal_test:
                    distance_traveled.append(aparticle.get_distance_traveled())
                    nhits = nhits + 1
        logger.info(f"nhits / nevents : {nhits} / {len(events)}")

        expected_mfp = (
            (self.nucleus.n_nucleons()+1)
            / (4.0/3.0*np.pi*self.nucleus.radius()**3.0)
            * settings().get('xsec')/10.0)**-1
        logger.info(f"Expected result: {expected_mfp} fm")

        _, ax = plt.subplots(1)
        _, scale = scipy.stats.expon.fit(distance_traveled)
        logger.info(f"Fitted result: {scale} fm")
        fit_label = rf"$\lambda$={scale:.2f} fm"
        sns.distplot(distance_traveled, ax=ax, kde=False, fit=scipy.stats.expon,
                     label='events', fit_kws={'label':fit_label})
        ax.set_xlabel(r'Distance traveled [fm]')
        ax.set_ylabel("Counts")
        ax.set_yscale('log')
        ax.set_title(
            f"Mean Free Path\n nhits / nevents : {nhits} / {len(events)}"
        )

        x_vals = np.linspace(0, 6, 1000)
        y_vals = np.exp(-x_vals/expected_mfp)/expected_mfp

        ax.plot(x_vals, y_vals, label='Expected')
        ax.legend()
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
        self.pid = 2212
        interaction = interactions.Interactions.create(settings().get_param('interaction'),
                                                       settings().get_param('interaction_file'))

        self.fsi = cascade.Cascade(interaction, self.cascade_prob)

    def generate_one_event(self):
        """ Generate one pN or nC event. """
        # Kick in a random direction
        # x in [0,1] --> cos(theta) in [-1,1] --> theta in [0,pi]
        # x in [0,1] --> phi in [0,2*pi]
        x = np.random.random(2)
        theta = np.arccos(2 * x[0] - 1)
        phi = 2 * np.pi * x[1]
        p_kick = settings().beam_energy
        self.nucleus.generate_config()
        particles = self.nucleus.nucleons()

        # Select a random particle to kick
        kicked_idx = np.random.choice(np.arange(len(particles)))
        self.fsi.set_kicked(kicked_idx)
        kicked_particle = particles[kicked_idx]
        kicked_particle.set_status(particle.ParticleStatus.internal_test)
        kicked_particle.set_momentum(vectors.Vector4(
            p_kick * np.sin(theta) * np.cos(phi),
            p_kick * np.sin(theta) * np.sin(phi),
            p_kick * np.cos(theta),
            np.sqrt(kicked_particle.mass()**2.0 + p_kick**2.0)))

        # Evolve the nucleus
        self.nucleus.set_nucleons(particles)
        self.fsi.mean_free_path(self.nucleus)

        return self.nucleus.nucleons()

    def finalize(self, events):
        """ Convert the events to a total cross-section. """
        nhits = 0
        for event in events:
            for aparticle in event:
                if aparticle.status() == particle.ParticleStatus.internal_test:
                    nhits = nhits + 1

        transparency = 1.0 - nhits / len(events)
        uncertainty = transparency / np.sqrt(nhits)
        logger.info(f"transparency: {transparency} +/- {uncertainty}")

        return transparency, uncertainty


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
        raise NotImplementedError

    def finalize(self, events):
        """ Convert the events to a total cross-section. """
        raise NotImplementedError
