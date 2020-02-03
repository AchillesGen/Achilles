""" Class to parse the nuchic input file and to store the settings. """

import shutil

import yaml
# from absl import logging

from .utils import make_path
#from .nucleus import Nucleus
from .histogram import Histogram


class _Settings:
    """ Class to read in the user settings from an input file or
    overwrite from commandline. """

    def __init__(self):
        self.nucleus = None

    def load(self, filename='run.yml'):
        """ Load a settings file. """
        try:
            with open(filename, 'r') as settings_file:
                self.__dict__.update(yaml.safe_load(settings_file))
        except FileNotFoundError:
            shutil.copyfile(make_path('template.yml'), 'run.yml')
            logging.fatal('{} could not be found. '
                          'Creating template at `run.yml` and '
                          'quitting.'.format(filename))

    @property
    def settings(self):
        """ Get all settings. """
        return self.__dict__

    def get(self, name, default=''):
        """ Get a given setting returning default if not found. """
        return _Settings.search(self.__dict__, name, default)[1]

    @staticmethod
    def search(dictionary, name, default=''):
        """ Search for a given setting. """
        found = False
        next_level = []
        for key, value in dictionary.items():
            if key == name:
                found = True
                return found, value
            if isinstance(value, dict):
                next_level.append(value)
        for sub_dict in next_level:
            found, result = _Settings.search(sub_dict, name, default)
            if found:
                return found, result
        return found, default


    # Run settings
    def get_run(self, name, default=''):
        """ Search for a setting in the run section. """
        return self.__dict__['run'].get(name, default)

    @property
    def run_settings(self):
        """ Return the dictionary of run settings. """
        return self.__dict__['run']

    @run_settings.setter
    def run_settings(self, name, value):
        """ Set a run setting. """
        self.__dict__['run'][name] = value

    @property
    def nevents(self):
        """ Return the requested number of generated events. """
        return self.run_settings['events']

    @property
    def beam_energy(self):
        """ Return the beam energy. """
        return self.run_settings['beam_energy']

    @property
    def angle(self):
        """ Return the angle of the outgoing electron in degrees. """
        return self.run_settings['angle']

    @property
    def cascade(self):
        """ Return if the cascade should be run. """
        return self.run_settings['cascade']

    @property
    def folding(self):
        """ Return if the folding function should be used. """
        return self.run_settings['folding']

    @property
    def output_format(self):
        """ Get the event output format. """
        return self.run_settings['output']

    # Nucleus settings
    def nucleus_params(self, name):
        return self.__dict__['{}_params'.format(name)]

    def nucleus_binding(self, name):
        return self.nucleus_params(name)['binding_energy']

    def nucleus_kf(self, name):
        return self.nucleus_params(name)['fermi_momentum']

    @property
    def nucleus_config(self):
        """ Get a dictionary for the nucleus configuration generation options. """
        return self.__dict__['nucleus_config']

    # Parameter settings
    def get_param(self, name, default=''):
        """ Search for a setting in the run section. """
        return self.__dict__['parameters'].get(name, default)

    @property
    def parameters(self):
        """ Return the dictionary of parameter settings. """
        return self.__dict__['parameters']

    @parameters.setter
    def parameters(self, name, value):
        """ Set a parameter value. """
        self.__dict__['parameters'][name] = value

    @property
    def distance(self):
        """ Maximum propagation distance of particles in cascade. """
        return self.parameters['cascade_distance']

    @property
    def folding_func(self):
        """ Get the user folding function. """
        return self.run_settings['folding_func']

    @property
    def config_type(self):
        """ Return the configuration type.

        Get the configuration type for how to setup the nucleus.
        Current options are either:
            - QMC: Quantum Monte Carlo configuration
            - MF: Mean field configuration
        """
        return self.parameters['config_type']

    # Other settings
    def get_histograms(self):
        """ Build the requested histograms from the yaml file. """
        histograms = {}
        for name, hist in self.__dict__['histograms'].items():
            histograms[name] = Histogram(**hist)

            # TODO: Store information on how to calculate

        return histograms


_SETTINGS = _Settings()


def settings():
    """ Accessor function for the settings. """
    return _SETTINGS
