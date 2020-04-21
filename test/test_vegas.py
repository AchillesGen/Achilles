""" Test vegas extension. """

import numpy as np

from nuchic.vegas import Integrator


def _func(x):
    return x[0] * x[1] * x[2]


def test_vegas_init():
    """ Test vegas sub-class initialization. """
    ndims = 10
    Integrator([[0, 1]]*ndims)


def test_vegas_weighted_events():
    """ Test vegas wgted events. """
    ndims = 3
    nevents = 1000
    vegas = Integrator([[0, 1]]*ndims)
    vegas_result = vegas(_func)

    points, weights, result = vegas.get_weighted_events(_func, nevents)

    assert len(points) == nevents
    assert len(weights) == nevents

    # Ensure answers are consistent at 2 sigma
    # (This should fail about 5% of the time)
    assert np.abs(vegas_result.mean - result.mean) \
        < 2*(vegas_result.sdev + result.sdev)


# def test_vegas_unweighted_events():
#     """ Test vegas wgted events. """
#     ndims = 3
#     nevents = 10000
#     vegas = Integrator([[0, 1]]*ndims)
#     result = vegas(_func)
#     v_result = result.mean
#     events = vegas.get_unweighted_events(_func, nevents)

#     assert len(events) == nevents

#     result = 0
#     result2 = 0
#     for event in events:
#         result += _func(event)
#         result2 += (_func(event))**2

#     assert np.abs(v_result - result) < np.sqrt(
#         (result2*nevents - result**2)/nevents)
