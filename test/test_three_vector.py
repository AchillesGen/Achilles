""" Tests for the three vector class. """

import pytest
import numpy as np
from nuchic.physics import Vector3 as Vec3


def test_vec3_init():
    """ Test Vec3 initialization. """
    vec1 = Vec3(1, 2, 3)
    vec2 = Vec3()
    vec2.set_xyz(1, 2, 3)
    vec3 = Vec3(3, 2, 1)

    assert vec1 == vec2
    assert vec1 != vec3


def test_vec3_getitem():
    """ Test Vec3 getitem. """
    vec = Vec3(1, 2, 3)
    assert vec[0] == 1
    assert vec[1] == 2
    assert vec[2] == 3
    with pytest.raises(Exception):
        print(vec[3])

    # assert isinstance(vec.array, np.ndarray)


def test_vec3_repr():
    """ Test Vec3 repr. """
    vec = Vec3(1, 2, 3)
    string = repr(vec)
    assert string == 'ThreeVector(1.000000, 2.000000, 3.000000)'


def test_vec3_str():
    """ Test Vec3 str. """
    vec = Vec3(1, 2, 3)
    string = str(vec)
    assert string == 'ThreeVector(1.000000, 2.000000, 3.000000)'


def test_vec3_add():
    """ Test Vec3 addition. """
    vec1 = Vec3(1, 2, 3)
    vec2 = Vec3(4, 3, 2)
    vec3 = Vec3(5, 5, 5)
    scalar = 3

    assert vec3 == vec1 + vec2
    with pytest.raises(Exception):
        print(vec1 + scalar)


def test_vec3_neg():
    """ Test Vec3 negation. """
    vec1 = Vec3(1, 2, 3)
    vec2 = Vec3(-1, -2, -3)

    assert vec1 == -vec2


def test_vec3_sub():
    """ Test Vec3 subtraction. """
    vec1 = Vec3(1, 2, 3)
    vec2 = Vec3(1, 2, 3)
    vec3 = Vec3()

    assert vec3 == vec1 - vec2


def test_vec3_mul():
    """ Test Vec3 multiplication. """
    scalar1 = 3
    vec1 = Vec3(1, 2, 3)
    vec2 = Vec3(3, 6, 9)
    scalar2 = 14

    assert vec1 * vec1 == scalar2
    assert vec1 * vec1 == vec1.dot(vec1)
    assert scalar1 * vec1 == vec1 * scalar1
    assert scalar1 * vec1 == vec2
    # assert vec1.__rmul__(vec2) == vec2 * vec1

    with pytest.raises(Exception):
        vec1.dot(scalar2)


def test_vec3_div():
    """ Test Vec3 division. """
    scalar = 3.0
    vec1 = Vec3(1, 2, 3)
    vec2 = Vec3(3, 6, 9)

    assert vec2 / scalar == vec1
    with pytest.raises(Exception):
        print(vec1 / vec2)


def test_vec3_magnitude():
    """ Test Vec3 magnitudes. """
    vec1 = Vec3(4, 3, 2)
    pvec2 = 29

    assert vec1.magnitude2() == pvec2
    assert vec1.magnitude() == np.sqrt(pvec2)


def test_vec3_angles():
    """ Test Vec3 angles. """
    vec = Vec3(0, 0, 4)
    vec1 = Vec3(1, 1, 1)

    assert vec.theta() == 0
    assert vec.phi() == 0
    assert vec1.phi() == np.pi / 4.0


def test_vec3_cross():
    """ Test Vec3 cross-product. """
    vec1 = Vec3(3, 2, 1)
    vec2 = Vec3(1, 2, 3)
    vec3 = Vec3(4, -8, 4)
    scalar = 5

    assert vec1.cross(vec2) == vec3
    assert vec1.cross(vec2) == -vec2.cross(vec1)

    with pytest.raises(Exception):
        vec1.cross(scalar)
