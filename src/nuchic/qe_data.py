""" Implement class to handle QE electron scattering data. """

import os
import pandas as pd

from .utils import make_path


class QuasielasticData:
    """ Class to hold and manipulate QE electron scattering data."""
    def __init__(self, element):
        data_file = make_path(element+'.dat', os.path.join('data', 'qe'))
        self.data = pd.read_csv(data_file,
                                sep=r'\s+',
                                names=('Z', 'A', 'Energy', 'Angle', 'Peak',
                                       'Data', 'StatUncertainty', 'citation'))
        self.element = element
        self.Z = self.data['Z'].values[0]
        self.A = self.data['A'].values[0]

    @property
    def citations(self):
        ''' Get all the experimental citations in the data set'''
        return self.data.citation.unique()

    @property
    def energies(self):
        ''' Get all the energy values in the data set'''
        return self.data.Energy.unique()

    @property
    def angles(self):
        ''' Get all the angle values in the data set'''
        return self.data.Angle.unique()

    @property
    def peaks(self):
        '''
            Get all the 4-momentum transfer at the top of the QE peak (x=1)
            values in the data set
        '''
        return self.data.Peak.unique()

    def get_data(self, energy, angle):
        """ Get the data for a specific energy and angle. """
        mask_energy = self.data['Energy'] == energy
        mask_angle = self.data['Angle'] == angle
        data_tmp = self.data[mask_energy & mask_angle]
        data_omega = data_tmp['Peak'].values
        data_dsigma = data_tmp['Data'].values
        data_error = data_tmp['StatUncertainty'].values

        return data_omega, data_dsigma, data_error
