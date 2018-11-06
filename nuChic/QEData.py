import pandas as pd

# Class to hold and manipulate QE electron scattering data
class QuasielasticData:
    def __init__(self, element):
        self.data = pd.from_csv(element+'.dat', sep='\s+', 
                names=('Z','A','Energy', 'Angle', 'Peak', 
                       'Data', 'StatUncertainty', 'citation'))
        self.element = element
        self.Z = self.data['Z'].values[0]
        self.A = self.data['A'].values[0]

    def GetCitations(self):
        ''' Get all the experimental citations in the data set'''
        return self.data.citation.unique()

    def GetEnergies(self):
        ''' Get all the energy values in the data set'''
        return self.data.Energy.unique()

    def GetAngles(self):
        ''' Get all the angle values in the data set'''
        return self.data.Angle.unique()

    def GetPeaks(self):
        ''' 
            Get all the 4-momentum transfer at the top of the QE peak (x=1) 
            values in the data set
        '''
        return self.data.Peak.unique()


