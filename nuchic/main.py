""" Main driving code for the nuChic program """

import os
from ast import literal_eval

from absl import flags, app, logging
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from .vegas import Integrator
from .inclusive import Quasielastic
from .nucleus import Nucleus
from .constants import MEV as MeV, FM
from .cascade import FSI
from .utils import momentum_sort
from .histogram import Histogram
from .input_parser import Settings

FLAGS = flags.FLAGS
flags.DEFINE_bool(
    'cascade', None, 'Flag to turn on/off the cascade', short_name='c')
flags.DEFINE_bool('folding', None, 'Flag to turn on/off the folding function')
flags.DEFINE_bool('testing', None, 'Flag to test the basics of the program')

DIR, FILE = os.path.split(__file__)


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
        self.settings = Settings(run_card)
        self.histograms = {'omega': Histogram([0, 500], 100),
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
                           'e_diff': Histogram([500, 1500], 100),
                           'escape': Histogram([-0.5, 12.5], 13),
                           'wgts': Histogram(bins=np.logspace(-4, 4, 400)),
                           }

        argon_nucleus = Nucleus(6, 12, 8.6*MeV, 225*MeV)
        self.inclusive = Quasielastic(0, self.settings,
                                      argon_nucleus, 730, 37.1)
        if FLAGS.cascade:
            self.fsi = FSI(argon_nucleus, self.settings.distance*FM)

        self.nevents = int(self.settings.nevents)

    def run(self):
        """ Runs the nuChic code with given settings."""
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

        if FLAGS.folding:
            wgt = wgt[1]

        self.fill_hists(momentum, variables, wgt*weight)

        return wgt*weight

    def fill_hists(self, momentum, variables, wgt):
        """ Fill the histograms. """
        if momentum is not None:
            self.histograms['omega'].fill(variables['omega'], wgt)
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
        # TODO: Use QuasielasticData class to load and plot the desired data
        # data = pd.read_csv(os.path.join(DIR, 'data', '12C.dat'), header=None,
        #                    sep=r'\s+',
        #                    names=[
        #                        'Z', 'A', 'Ee', 'Angle', 'omega',
        #                        'dsigma', 'error', 'citation'
        #                    ])

        # energy = 0.730
        # angle = 37.1

        # mask_energy = data['Ee'] == energy
        # mask_angle = data['Angle'] == angle
        # data_tmp = data[mask_energy & mask_angle]
        # data_omega = data_tmp['omega'].values
        # data_dsigma = data_tmp['dsigma'].values
        # data_error = data_tmp['error'].values

        # bins = np.linspace(0, data_omega[-2]*1000, data_omega[-2]*1000/5.+1)
        # bin_widths = np.diff(bins)
        # print(bins)

        # hist, bins = np.histogram(omega, weights=wgts, bins=bins)
        # hist_f, bins = np.histogram(omega, weights=wgts_f, bins=bins)

        # hist /= bin_widths
        # hist_f /= bin_widths

        _, ax1 = plt.subplots(nrows=1, ncols=1)
        self.histograms['omega'].plot(ax1)
        # ax1.errorbar(data_omega*1000, data_dsigma,
        #              yerr=data_error, fmt='o', color='red')
        # # ax1.plot(omega_f,dsigma, color='red', ls='steps')
        # ax1.plot(bins[:-1], hist, color='blue', ds='steps', label='No FSI')
        # ax1.plot(bins[:-1], hist_f, color='green', ds='steps', label='FSI')
        # ax1.legend()
        # ax1.hist(omega,weights=wgts,bins=data_omega[:-1]*1000,color='blue')
        # ax1.hist(omega,weights=wgts_f,bins=data_omega[:-1]*1000,color='green')
        # ax2.hist(mom,weights=wgts,bins=400)
        # ax3.hist(energies,weights=wgts,bins=400)

        fig2, axs = plt.subplots(nrows=2, ncols=2)
        self.histograms['e_post'].plot(axs[0][0])
        self.histograms['px_post'].plot(axs[0][1])
        self.histograms['py_post'].plot(axs[1][0])
        self.histograms['pz_post'].plot(axs[1][1])
        fig2.suptitle('After Cascade')

        fig3, axs2 = plt.subplots(nrows=2, ncols=2)
        self.histograms['e_pre'].plot(axs2[0][0])
        self.histograms['px_pre'].plot(axs2[0][1])
        self.histograms['py_pre'].plot(axs2[1][0])
        self.histograms['pz_pre'].plot(axs2[1][1])
        fig3.suptitle('Before Cascade')

        fig4, axs3 = plt.subplots(nrows=2, ncols=2)
        self.histograms['e_diff'].plot(axs3[0][0])
        self.histograms['px_diff'].plot(axs3[0][1])
        self.histograms['py_diff'].plot(axs3[1][0])
        self.histograms['pz_diff'].plot(axs3[1][1])
        fig4.suptitle('Difference in momentum')

        if FLAGS.cascade:
            _, ax8 = plt.subplots(nrows=1, ncols=1)
            self.histograms['escape'].plot(ax8)
            ax8.set_yscale('log')

        plt.show()


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
