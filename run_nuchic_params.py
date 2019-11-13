""" Run nuchic over a range of parameters. """

import os
import subprocess
import shutil
from string import Template

import numpy as np

from nuchic.qe_data import QuasielasticData

PARAMETERS = ['beam_energy', 'angle']


def get_data():
    qe_data = QuasielasticData(str('12C'))
    values = qe_data.data.groupby(['Energy', 'Angle']).groups.keys()
    return values


def create_run(parameters):
    """ Create a run card from a dictionary of parameters. """
    template = open('run_template.yml')
    src = Template(template.read())
    template.close()
    result = src.substitute(parameters)
    with open('run.yml', 'w') as run_card:
        for line in result:
            run_card.write(line)


def main():
    """ Loop over all the parameters and generate results. """
    if not os.path.exists('Results'):
        os.mkdir('Results')

    values = get_data()

    for value in values:
        energy = int(np.rint(value[0]*1000))
        angle = np.rint(value[1]*100)/100
        print(energy, angle)
        create_run({PARAMETERS[0]: energy, PARAMETERS[1]: angle})

        subprocess.call('nuchic')

        shutil.move('omega.png',
                    os.path.join('Results',
                                 'omega_{}_{}.png'.format(energy, angle)))

        shutil.move('before_cascade.png',
                    os.path.join('Results',
                                 'before_cascade_{}_{}.png'.format(energy,
                                                                   angle)))


if __name__ == '__main__':
    main()
