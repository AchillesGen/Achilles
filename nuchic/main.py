""" Main driving code for the nuChic program """

import os
from ast import literal_eval

from absl import flags, app, logging
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from tqdm import tqdm

from .vegas import Integrator
from .inclusive import Quasielastic
from .nucleus import Nucleus
from .constants import MEV as MeV
from .cascade import FSI
from .utils import momentum_sort
from .histogram import Histogram

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
    def __init__(self, nevents):
        self.histograms = {'omega': Histogram([0, 500], 500),
                           'px_pre': Histogram([-500, 500], 1000),
                           'py_pre': Histogram([-500, 500], 1000),
                           'pz_pre': Histogram([-500, 500], 1000),
                           'e_pre': Histogram([-500, 500], 1000),
                           'px_post': Histogram([-500, 500], 1000),
                           'py_post': Histogram([-500, 500], 1000),
                           'pz_post': Histogram([-500, 500], 1000),
                           'e_post': Histogram([-500, 500], 1000),
                           'px_diff': Histogram([-500, 500], 1000),
                           'py_diff': Histogram([-500, 500], 1000),
                           'pz_diff': Histogram([-500, 500], 1000),
                           'e_diff': Histogram([-500, 500], 1000),
                           'escape': Histogram([-0.5, 12.5], 13),
                           'wgts': Histogram(bins=np.logspace(-4, 4, 400)),
                           }

        argon_nucleus = Nucleus(6, 12, 8.6*MeV, 225*MeV)
        self.inclusive = Quasielastic(argon_nucleus, 730, 37.1)
        if FLAGS.cascade:
            self.fsi = FSI(argon_nucleus, 1)

        self.nevents = int(nevents)

    def run(self):
        """ Runs the nuChic code with given settings."""
        ndims = 4
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
        # plot_results(self.histograms)

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

                    self.histograms['e_pre'].fill(momentum_post.energy, wgt)
                    self.histograms['px_pre'].fill(momentum_post.p_x, wgt)
                    self.histograms['py_pre'].fill(momentum_post.p_y, wgt)
                    self.histograms['pz_pre'].fill(momentum_post.p_z, wgt)

                    self.histograms['e_diff'].fill(
                        momentum_post.energy - momentum.energy, wgt)
                    self.histograms['px_diff'].fill(
                        momentum_post.p_x - momentum.p_x, wgt)
                    self.histograms['py_diff'].fill(
                        momentum_post.p_y - momentum.p_y, wgt)
                    self.histograms['pz_diff'].fill(
                        momentum_post.p_z - momentum.p_z, wgt)


def plot_results(histograms):
    """ Plot results of the run. """
    data = pd.read_csv(os.path.join(DIR, 'data', '12C.dat'), header=None,
                       sep=r'\s+',
                       names=[
                           'Z', 'A', 'Ee', 'Angle', 'omega',
                           'dsigma', 'error', 'citation'
                       ])

    energy = 0.730
    angle = 37.1

    mask_energy = data['Ee'] == energy
    mask_angle = data['Angle'] == angle
    data_tmp = data[mask_energy & mask_angle]
    data_omega = data_tmp['omega'].values
    data_dsigma = data_tmp['dsigma'].values
    data_error = data_tmp['error'].values

    bins = np.linspace(0, data_omega[-2]*1000, data_omega[-2]*1000/5.+1)
    bin_widths = np.diff(bins)
    print(bins)

    # hist, bins = np.histogram(omega, weights=wgts, bins=bins)
    # hist_f, bins = np.histogram(omega, weights=wgts_f, bins=bins)

    # hist /= bin_widths
    # hist_f /= bin_widths

    # _, ax1 = plt.subplots(nrows=1, ncols=1)
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
    axs[0][0].hist(p_e, weights=wgts, bins=400)
    axs[0][1].hist(p_px, weights=wgts, bins=400)
    axs[1][0].hist(p_py, weights=wgts, bins=400)
    axs[1][1].hist(p_pz, weights=wgts, bins=400)
    fig2.suptitle('After Cascade')

    fig3, axs2 = plt.subplots(nrows=2, ncols=2)
    axs2[0][0].hist(p_e_pre, weights=wgts, bins=400)
    axs2[0][1].hist(p_px_pre, weights=wgts, bins=400)
    axs2[1][0].hist(p_py_pre, weights=wgts, bins=400)
    axs2[1][1].hist(p_pz_pre, weights=wgts, bins=400)
    fig3.suptitle('Before Cascade')

    fig4, axs3 = plt.subplots(nrows=2, ncols=2)
    axs3[0][0].hist(p_e_diff, weights=wgts, bins=400)
    axs3[0][1].hist(p_px_diff, weights=wgts, bins=400)
    axs3[1][0].hist(p_py_diff, weights=wgts, bins=400)
    axs3[1][1].hist(p_pz_diff, weights=wgts, bins=400)
    fig4.suptitle('Difference in momentum')

#    fig5, ax7 = plt.subplots(nrows=1, ncols=1)
#    ax7.hist(wgts,bins=np.logspace(np.log10(np.min(wgts)),
#                                   np.log10(np.max(wgts)), 100))
#    ax7.set_xscale('log')
#    ax7.set_yscale('log')
#
#    print(np.mean(wgts)/np.max(wgts))

    if FLAGS.cascade:
        _, ax8 = plt.subplots(nrows=1, ncols=1)
        ax8.hist(escape, weights=wgts, bins=np.linspace(-0.5, 10.5, 12))
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
        nuchic = NuChic(1e5)
        nuchic.run()

    logging.info('nuChic finished successfully!')


def nu_chic():
    """ Entry point for running the program """
    app.run(main)


if __name__ == "__main__":
    app.run(main)
