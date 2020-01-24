""" Main driver code. """

from tqdm import tqdm

import nuchic.run_modes as run_modes
from .config import settings


class NuChic:
    """ Main driver class for the nuchic code. """
    def __init__(self, run_card):
        settings().load(run_card)
        self.calc = run_modes.RunMode(settings().get_run('mode', 'generate'), 'test')
        self.nevents = settings().nevents

    def run(self):
        """ Generate events for the given input settings. """
        for _ in tqdm(range(self.nevents), ncols=80):
            event = self.calc.generate_one_event()


def nuchic():
    """ External entry point. """
    NuChic('run.yml')


if __name__ == '__main__':
    nuchic()
