""" Implement tests of nuChic utilities. """

import numpy as np

import nuchic.utils as utils


def test_transforms():
    """ Test coordinate transformations. """
    coords = np.random.random((3, 20))
    coords[0] *= 10
    coords[1] *= np.pi
    coords[2] *= 2*np.pi

    x_coords = utils.to_cartesian(coords)
    r_coords = utils.to_spherical(x_coords)
    assert x_coords.shape == coords.shape
    assert np.allclose(r_coords, coords)

    coords2 = coords.T
    x_coords2 = utils.to_cartesian(coords2, axis=1)
    r_coords2 = utils.to_spherical(x_coords2, axis=1)
    assert x_coords2.shape == coords2.shape
    assert np.allclose(x_coords2.T, x_coords)
    assert np.allclose(r_coords2, coords2)
