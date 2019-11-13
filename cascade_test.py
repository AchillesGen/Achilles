""" Run nuchic over a range of parameters. """

import os
import subprocess
from string import Template

import numpy as np

PARAMETERS = ['beam_energy']


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

    emin = 250
    emax = 1000
    energies = np.linspace(emin, emax, (emax-emin)/10+1)

    for energy in energies:
        print(energy)
        create_run({PARAMETERS[0]: energy})
        subprocess.call(['nuchic', '--cascade_test'])


if __name__ == '__main__':
    main()
