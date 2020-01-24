""" Implementation of utility functions for the nuchic code.

Some useful utilities / decorators for debuging
"""

import os
import time
import functools
import logging
import numpy as np
from scipy.special import gammainc

LOGGER = logging.getLogger(__name__)
DIR, FILE = os.path.split(__file__)


def timing(fcn):
    """Time the execution of fcn. Use as decorator."""
    @functools.wraps(fcn)
    def wrap(*args, **kwargs):
        """Wrapped version of the function."""
        t_initial = time.perf_counter()
        result = fcn(*args, **kwargs)
        t_final = time.perf_counter()
        LOGGER.debug(
            "TIMING: %s took: %.4e sec",
            fcn.__name__,
            t_final - t_initial
        )
        return result
    return wrap


def make_path(filename, path=''):
    """ Make a file path to an installed nuchic file. """
    loc = os.path.join(DIR, path, filename)
    return loc


def momentum_sort(elem):
    """ Sort the array by momentum. """
    return elem.mom.mom


@timing
def to_cartesian(coords, axis=0):
    """Convert spherical coordinates to cartesian coordinates. """
    coords = np.moveaxis(coords, axis, 0)
    r, theta, phi = coords
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    coords = np.moveaxis(coords, 0, axis)
    return np.moveaxis(np.array([x, y, z]), 0, axis)


def to_spherical(coords, axis=0):
    """Convert cartesian coordinates to spherical coordinates. """
    coords = np.moveaxis(coords, axis, 0)
    x, y, z = coords
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    phi = np.where(phi > 0, phi, phi + 2*np.pi)
    coords = np.moveaxis(coords, 0, axis)
    return np.moveaxis(np.array([r, theta, phi]), 0, axis)

def rand_sphere(radius, num):
    """ Generate points uniformly within a sphere of a given radius.

    Args:
        radius (float): Radius of the sphere to generate within
        num (int): Number of points to generate within the sphere

    Returns:
        nd.ndarray: A (num, 3) sized array of the randomly sampled points
    """
    pts = np.random.normal(size=(num, 3))
    ssq = np.sum(pts**2, axis=1)
    tmp = radius*gammainc(3./2., ssq/2.)**(1./3.)/np.sqrt(ssq)
    tmp_tiled = np.tile(tmp.reshape(num, 1), (1, 3))
    return np.multiply(pts, tmp_tiled)
