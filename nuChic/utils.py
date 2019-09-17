""" Implementation of utility functions for the nuChic code. """

import numpy as np


def momentum_sort(elem):
    """ Sort the array by momentum. """
    return elem.mom.P()


def to_cartesian(coords):
    """Convert spherical coordinates to cartesian coordinates. """
    # r, theta, phi = coords
    r = coords[:, 0]
    theta = coords[:, 1]
    phi = coords[:, 2]
    x = r * np.sin(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.cos(phi)
    z = r * np.cos(theta)
    return np.transpose(np.array([x, y, z]))
