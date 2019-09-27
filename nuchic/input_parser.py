""" Class to parse the nuchic input file and to store the settings. """

import yaml

from nuchic.histogram import Histogram


class Settings:
    """ Class to read in the user settings from an input file or
    overwrite from commandline. """
    def __init__(self, filename='run.yml'):
        with open(filename, 'r') as settings_file:
            self.settings = yaml.safe_load(settings_file)

        # TODO: Override settings with commandline arguments if given

    @property
    def nevents(self):
        """ Return the requested number of generated events. """
        return self.settings['run']['events']

    @property
    def cascade(self):
        """ Return if the cascade should be run. """
        return self.settings['run']['cascade']

    @property
    def folding(self):
        """ Return if the folding function should be used. """
        return self.settings['run']['folding']

    @property
    def output_format(self):
        """ Get the event output format. """
        return self.settings['run']['output']

    @property
    def distance(self):
        """ Maximum propagation distance of particles in cascade. """
        return self.settings['parameters']['cascade_distance']

    def get_histograms(self):
        """ Build the requested histograms from the yaml file. """
        histograms = {}
        for name, hist in self.settings['histograms'].items():
            histograms[name] = Histogram(**hist)

            # TODO: Store information on how to calculate

        return histograms
