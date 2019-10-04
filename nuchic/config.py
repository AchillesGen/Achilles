""" Class to parse the nuchic input file and to store the settings. """

import yaml

from nuchic.histogram import Histogram


class _Settings:
    """ Class to read in the user settings from an input file or
    overwrite from commandline. """

    def __init__(self):
        pass

    def load(self, filename='run.yml'):
        """ Load a settings file. """
        with open(filename, 'r') as settings_file:
            self.__dict__.update(yaml.safe_load(settings_file))

#    def __getattr__(self, name):
#        return self.__dict__.get(name, False)

    @property
    def settings(self):
        """ Get all settings. """
        return self.__dict__

    def search(self, name):
        """ Search for a given setting. """

    @property
    def run_settings(self):
        """ Return the dictionary of run settings. """
        return self.__dict__['run']

    @run_settings.setter
    def run_settings(self, name, value):
        """ Set a run setting. """
        self.__dict__['run'][name] = value

    @property
    def parameters(self):
        """ Return the dictionary of parameter settings. """
        return self.__dict__['parameters']

    @parameters.setter
    def parameters(self, name, value):
        """ Set a parameter value. """
        self.__dict__['parameters'][name] = value

    @property
    def nevents(self):
        """ Return the requested number of generated events. """
        return self.run_settings['events']

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

    @property
    def distance(self):
        """ Maximum propagation distance of particles in cascade. """
        return self.parameters['cascade_distance']

    @property
    def folding_func(self):
        """ Get the user folding function. """
        return self.run_settings['folding_func']

    def get_histograms(self):
        """ Build the requested histograms from the yaml file. """
        histograms = {}
        for name, hist in self.__dict__['histograms'].items():
            histograms[name] = Histogram(**hist)

            # TODO: Store information on how to calculate

        return histograms


SETTINGS = _Settings()
