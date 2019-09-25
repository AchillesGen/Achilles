""" Implement class to handle QE electron scattering data. """

import pandas as pd


class QuasielasticData:
    """ Class to hold and manipulate QE electron scattering data."""
    def __init__(self, element):
        self.data = pd.from_csv(element+'.dat', sep=r'\s+',
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
