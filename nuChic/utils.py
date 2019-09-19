""" Implementation of utility functions for the nuChic code. """

import numpy as np


def momentum_sort(elem):
    """ Sort the array by momentum. """
    return elem.mom.mom


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
