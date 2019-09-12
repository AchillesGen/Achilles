""" Main driving code for the nuChic program """

import os

from absl import flags, app, logging
import vegas
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from .inclusive import Quasielastic
from .Nucleus import Nucleus
from .Constants import MeV
from .cascade import FSI
from .utils import momentum_sort

FLAGS = flags.FLAGS
flags.DEFINE_bool(
    'cascade', None, 'Flag to turn on/off the cascade', short_name='c')
flags.DEFINE_bool('folding', None, 'Flag to turn on/off the folding function')
flags.DEFINE_bool('testing', None, 'Flag to test the basics of the program')

DIR, FILE = os.path.split(__file__)


def run():
    """ Runs the nuChic code with given settings."""
    omega = []

    p_px_pre = []
    p_py_pre = []
    p_pz_pre = []
    p_e_pre = []

    p_px = []
    p_py = []
    p_pz = []
    p_e = []

    p_px_diff = []
    p_py_diff = []
    p_pz_diff = []
    p_e_diff = []

    wgts = []
    wgts_f = []
    escape = []

    ndims = 3
    if FLAGS.folding:
        ndims += 1

    argon_nucleus = Nucleus(6, 12, 8.6*MeV, 225*MeV)
    inclusive = Quasielastic(argon_nucleus, 730, 37.1)
    if FLAGS.cascade:
        fsi = FSI(argon_nucleus, 1)

    integ = vegas.Integrator([[0, 1]]*ndims)

    result = integ(inclusive.GenerateWeight, nitn=20, neval=1e4)

    print(result.summary())

    nitn = 20
    count = 0
    zeros = 0
    zeros_f = 0
    for _ in range(nitn):
        for point, weight in integ.random():

            if count % 10000 == 0:
                print(count)
            count += 1

            wgt = inclusive.GenerateWeight(point)

            if FLAGS.folding:
                wgt_f = wgt[0]
                wgt = wgt[1]
            else:
                wgt_f = 0

            if wgt == 0:
                zeros += 1

            if wgt_f == 0:
                zeros_f += 1

            # mom.append(p_int)
            # energies.append(e_int)
            momentum = inclusive.GenerateMomentum()
            if momentum is not None:
                omega.append(inclusive.w)
                p_e_pre.append(momentum.E)
                p_px_pre.append(momentum.px)
                p_py_pre.append(momentum.py)
                p_pz_pre.append(momentum.pz)
                wgts.append(wgt*weight/nitn)
                wgts_f.append(wgt_f*weight/nitn)

                if FLAGS.cascade:
                    fsi.kick(momentum)
                    escaped_part = fsi()
                    escape.append(len(escaped_part))
                    fsi.reset()

                    if not escaped_part:
                        escaped_part.sort(reverse=True, key=momentum_sort)
                        momentum = escaped_part[0].mom

                p_e.append(momentum.E)
                p_px.append(momentum.px)
                p_py.append(momentum.py)
                p_pz.append(momentum.pz)

                p_e_diff.append(momentum.E-p_e_pre[-1])
                p_px_diff.append(momentum.px-p_px_pre[-1])
                p_py_diff.append(momentum.py-p_py_pre[-1])
                p_pz_diff.append(momentum.pz-p_pz_pre[-1])

    print(sum(wgts), sum(wgts_f))
    print(count)
    print(zeros)
    print(zeros_f)

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

    hist, bins = np.histogram(omega, weights=wgts, bins=bins)
    hist_f, bins = np.histogram(omega, weights=wgts_f, bins=bins)

    hist /= bin_widths
    hist_f /= bin_widths

    _, ax1 = plt.subplots(nrows=1, ncols=1)
    ax1.errorbar(data_omega*1000, data_dsigma,
                 yerr=data_error, fmt='o', color='red')
    # ax1.plot(omega_f,dsigma, color='red', ls='steps')
    ax1.plot(bins[:-1], hist, color='blue', ds='steps', label='No FSI')
    ax1.plot(bins[:-1], hist_f, color='green', ds='steps', label='FSI')
    ax1.legend()
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
        run()

    logging.info('nuChic finished successfully!')


def nu_chic():
    """ Entry point for running the program """
    app.run(main)


if __name__ == "__main__":
    app.run(main)
