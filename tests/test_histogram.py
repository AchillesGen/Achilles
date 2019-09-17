""" Testing histogram code. """

import numpy as np
import pytest

from nuChic.histogram import Histogram


def test_histogram_init():
    """ Test histogram initialization. """
    bins = np.linspace(0., 1., 101)
    hist = Histogram(bins=bins)
    hist2 = Histogram(min_val=0., max_val=1., bins=101)

    assert np.allclose(hist.bins, hist2.bins)
    assert hist.min == hist2.min
    assert hist.max == hist2.max

    bins = np.logspace(np.log10(1e-3), np.log10(1.), 101)
    hist = Histogram(bins=bins)
    hist2 = Histogram(min_val=1e-3, max_val=1., bins=101, scale='log')

    assert np.allclose(hist.bins, hist2.bins)
    assert hist.min == hist2.min
    assert hist.max == hist2.max

    with pytest.raises(TypeError):
        print(Histogram(min_val=0., max_val=1., bins=bins))

    with pytest.raises(ValueError):
        print(Histogram(min_val=0., max_val=1., bins=101, scale='fjskla'))


def test_histogram_fill():
    """ Test filling of histogram. """
    bins = np.linspace(0., 100., 101)
    hist = Histogram(bins=bins)
    hist.fill(0.5, 10)
    hist.fill(1.5, 2)
    hist.fill(85.5, 1405)

    assert hist.weights[np.searchsorted(bins, 0.5)-1] == 10
    assert hist.weights[np.searchsorted(bins, 1.5)-1] == 2
    assert hist.weights[np.searchsorted(bins, 85.5)-1] == 1405

    hist.fill(-5, 10)
    hist.fill(200, 10)

    assert hist.underflow == 10
    assert hist.overflow == 10

    wgts = np.array(hist.weights)
    hist.fill(2341)
    assert np.allclose(wgts, hist.weights)


def test_histogram_str():
    """ Test string output of histogram. """
    hist = Histogram(0., 1., 2)
    line_format = '{:^15}{:^15}{:^15}{:^15}\n'
    string = line_format.format('left', 'right', 'wgt', 'wgt2')
    string += line_format.format(0.0, 1.0, 0.0, 0.0)

    assert str(hist) == string


def test_histogram_integral():
    """ Test the integration of the histogram. """
    hist = Histogram(0., 1., 3)
    hist.fill(0.25, 10)
    hist.fill(0.75, 20)

    assert hist.integral() == 0.5*10+0.5*20
