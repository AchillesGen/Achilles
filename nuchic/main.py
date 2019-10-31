""" Main driving code for the nuChic program """
from ast import literal_eval
import os

from absl import flags, app, logging
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from .vegas import Integrator
from .inclusive import Quasielastic
from .cascade import FSI

from .constants import FM, GEV, MQE
from .utils import momentum_sort
from .histogram import Histogram
from .config import settings
from .qe_data import QuasielasticData
from .particle import Particle
from .four_vector import Vec4
from .three_vector import Vec3

FLAGS = flags.FLAGS
flags.DEFINE_bool(
    'cascade', None, 'Flag to turn on/off the cascade', short_name='c')
flags.DEFINE_bool('folding', None, 'Flag to turn on/off the folding function')
flags.DEFINE_bool('cascade_test', None, 'Flag to test the cascade code')
flags.DEFINE_bool('testing', None, 'Flag to test the basics of the program')


# TODO: Figure out how to implement this correctly
def init_histograms(filename='hists.txt'):
    """ Load the desired histograms to be plotted from input file. """
    histograms = {}
    with open(filename) as hist_file:
        for line in hist_file:
            print(line)
            tmp = literal_eval(line)
            histograms[tmp.name] = tmp


class NuChic:
    """ Main driver class for the nuChic code. """
    def __init__(self, run_card):
        settings().load(run_card)
        self.histograms = {'omega': Histogram([0, 0.6], 200),
                           'omega2': Histogram([0, 0.6], 200),
                           'px_pre': Histogram([-500, 500], 100),
                           'py_pre': Histogram([-500, 500], 100),
                           'pz_pre': Histogram([-500, 500], 100),
                           'e_pre': Histogram([500, 1500], 100),
                           'px_post': Histogram([-500, 500], 100),
                           'py_post': Histogram([-500, 500], 100),
                           'pz_post': Histogram([-500, 500], 100),
                           'e_post': Histogram([500, 1500], 100),
                           'px_diff': Histogram([-500, 500], 100),
                           'py_diff': Histogram([-500, 500], 100),
                           'pz_diff': Histogram([-500, 500], 100),
                           'e_diff': Histogram([-500, 500], 100),
                           'escape': Histogram([-0.5, 12.5], 13),
                           'wgts': Histogram(bins=np.logspace(-4, 4, 400)),
                           }

        self.inclusive = Quasielastic(0,
                                      settings().beam_energy,
                                      settings().angle)
        if FLAGS.cascade:
            self.fsi = FSI(settings().distance*FM)

        self.nevents = int(settings().nevents)

    def run(self):
        """ Runs the nuChic code with given settings."""

        if FLAGS.cascade_test:
            print(self.calc_cross_section())
            return

        ndims = 3
        if FLAGS.folding:
            ndims += 1

        integ = Integrator([[0, 1]]*ndims)
        result = integ(self.generate_one_event, nitn=10, neval=1e4)
        print(result.summary())

        zeros = 0

        points, weights, _ = integ.generate(self.nevents)
        for i in tqdm(range(self.nevents), ncols=80):
            if self.generate_one_event(points[i], weights[i]) == 0:
                zeros += 1

        print(zeros)
        self.plot_results()

    def generate_one_event(self, point, weight=None):
        """ Generate the output for a single event. """

        wgt, variables, qval = self.inclusive.generate_weight(point)
        if weight is None:
            return wgt

        if wgt == 0:
            return 0

        momentum = self.inclusive.generate_momentum(variables, qval)

        self.fill_hists(momentum, variables, wgt*weight)

        return wgt*weight

    def fill_hists(self, momentum, variables, wgt):
        """ Fill the histograms. """
        if momentum is not None:
            self.histograms['omega'].fill(variables.omega / GEV, wgt)
            self.histograms['e_pre'].fill(momentum.energy, wgt)
            self.histograms['px_pre'].fill(momentum.p_x, wgt)
            self.histograms['py_pre'].fill(momentum.p_y, wgt)
            self.histograms['pz_pre'].fill(momentum.p_z, wgt)
            self.histograms['wgts'].fill(wgt)

            if FLAGS.cascade:
                self.fsi.kick(momentum)
                escaped_part = self.fsi()
                self.histograms['escape'].fill(len(escaped_part), wgt)
                self.fsi.reset()

                if escaped_part:
                    self.histograms['omega2'].fill(variables.omega / GEV, wgt)
                    escaped_part.sort(reverse=True, key=momentum_sort)
                    momentum_post = escaped_part[0].mom

                    self.histograms['e_post'].fill(momentum_post.energy, wgt)
                    self.histograms['px_post'].fill(momentum_post.p_x, wgt)
                    self.histograms['py_post'].fill(momentum_post.p_y, wgt)
                    self.histograms['pz_post'].fill(momentum_post.p_z, wgt)

                    self.histograms['e_diff'].fill(
                        momentum_post.energy - momentum.energy, wgt)
                    self.histograms['px_diff'].fill(
                        momentum_post.p_x - momentum.p_x, wgt)
                    self.histograms['py_diff'].fill(
                        momentum_post.p_y - momentum.p_y, wgt)
                    self.histograms['pz_diff'].fill(
                        momentum_post.p_z - momentum.p_z, wgt)

    def plot_results(self):
        """ Plot results of the run. """
        qe_data = QuasielasticData(str(settings().nucleus))
        energy = settings().beam_energy / GEV
        angle = settings().angle

        data = qe_data.get_data(energy, angle)

        _, ax1 = plt.subplots(nrows=1, ncols=1)
        self.histograms['omega'].plot(ax1)
        self.histograms['omega2'].plot(ax1)
        ax1.errorbar(data[0], data[1],
                     yerr=data[2], fmt=',', color='red',
                     capsize=2.5)

        plt.savefig('omega.png', bbox_inches='tight')

        fig3, axs2 = plt.subplots(nrows=2, ncols=2)
        self.histograms['e_pre'].plot(axs2[0][0])
        self.histograms['px_pre'].plot(axs2[0][1])
        self.histograms['py_pre'].plot(axs2[1][0])
        self.histograms['pz_pre'].plot(axs2[1][1])
        fig3.suptitle('Before Cascade')

        plt.savefig('before_cascade.png', bbox_inches='tight')

        if FLAGS.cascade:
            fig2, axs = plt.subplots(nrows=2, ncols=2)
            self.histograms['e_post'].plot(axs[0][0])
            self.histograms['px_post'].plot(axs[0][1])
            self.histograms['py_post'].plot(axs[1][0])
            self.histograms['pz_post'].plot(axs[1][1])
            fig2.suptitle('After Cascade')

            plt.savefig('after_cascade.png', bbox_inches='tight')

            fig4, axs3 = plt.subplots(nrows=2, ncols=2)
            self.histograms['e_diff'].plot(axs3[0][0])
            self.histograms['px_diff'].plot(axs3[0][1])
            self.histograms['py_diff'].plot(axs3[1][0])
            self.histograms['pz_diff'].plot(axs3[1][1])
            fig4.suptitle('Difference in momentum')

            plt.savefig('difference_cascade.png', bbox_inches='tight')

            _, ax8 = plt.subplots(nrows=1, ncols=1)
            self.histograms['escape'].plot(ax8)
            ax8.set_yscale('log')

            plt.savefig('escape.png', bbox_inches='tight')

    def calc_cross_section(self):
        """ Calculate the pC and nC cross-sections

        Returns: pC and nC cross-sections

        """
        # Initialize FSI
        self.fsi = FSI(settings().distance*FM, final=False)

        # Initialize variables
        radius = settings().nucleus.radius
        xsec = np.pi*radius**2
        nscatter_p = 0
        nscatter_n = 0

        # Loop over events
        for _ in range(self.nevents):
            position_r = np.random.rand(1)[0]*radius
            position_th = np.random.rand(1)[0]*2*np.pi
            position = Vec3(position_r*np.cos(position_th),
                            position_r*np.sin(position_th),
                            -3*radius)
            momentum = Vec4(np.sqrt(MQE**2 + settings().beam_energy**2),
                            0, 0,
                            settings().beam_energy)

            # Perform proton calculation
            proton = Particle(2212, momentum, position)
            self.fsi.kicked_idxs = [len(self.fsi.nucleons)]
            self.fsi.nucleons.append(proton)
            self.fsi.nucleons[-1].status = -1
            self.fsi()
            if self.fsi.scatter:
                nscatter_p += 1
            self.fsi.reset()

            # Perform neutron calculation
            neutron = Particle(2112, momentum, position)
            self.fsi.kicked_idxs = [len(self.fsi.nucleons)]
            self.fsi.nucleons.append(neutron)
            self.fsi.nucleons[-1].status = -1
            self.fsi()
            if self.fsi.scatter:
                nscatter_n += 1
            self.fsi.reset()

        return xsec * nscatter_p/self.nevents, xsec * nscatter_n/self.nevents


def main(argv):
    '''
    The Main prgram code
    Authors:
    '''
    del argv

    logging.info('Starting nuChic...')
    if not FLAGS.testing:
        nuchic = NuChic('run.yml')
        nuchic.run()

    logging.info('nuChic finished successfully!')


def nu_chic():
    """ Entry point for running the program """
    app.run(main)


if __name__ == "__main__":
    app.run(main)
