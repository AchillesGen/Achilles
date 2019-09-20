""" Test vegas extension. """

import numpy as np

from nuChic.vegas import Integrator


def _func(x):
    return x[0] * x[1] * x[2]


def test_vegas_init():
    """ Test vegas sub-class initialization. """
    ndims = 10
    Integrator([[0, 1]]*ndims)


def test_vegas_events():
    """ Test vegas wgted events. """
    ndims = 3
    nevents = 10000
    vegas = Integrator([[0, 1]]*ndims)
    result = vegas(_func, neval=1e6)
    print(result.summary())
    vegas.set(beta=-1)
    points, weights = vegas.get_events(nevents)

    assert len(points) == nevents
    assert len(weights) == nevents

    result = 0
    for i, _ in enumerate(points):
        result += _func(points[i]) * weights[i]

    print(result)

    assert False
