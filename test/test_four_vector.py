""" Tests Four Vector."""
import numpy as np
import pytest

from nuchic.physics import Vector4


def test_vec4_init():
    """ Test four vector initialization. """
    vec1 = Vector4(2, 3, 4, 1)
    vec2 = Vector4()
    vec2.set_pxpypze(2, 3, 4, 1)
    vec3 = Vector4(4, 3, 2, 1)

    assert vec1 == vec2
    assert vec1 != vec3


def test_vec4_getitem():
    """ Test four vector getitem. """
    vec = Vector4(1, 2, 3, 4)
    assert vec[0] == 1
    assert vec[1] == 2
    assert vec[2] == 3
    assert vec[3] == 4
    with pytest.raises(Exception):
        print(vec[4])


def test_vec4_repr():
    """ Test repr output. """
    vec = Vector4(1, 2, 3, 4)
    string = repr(vec)
    assert string == 'FourVector(1.000000, 2.000000, 3.000000, 4.000000)'


def test_vec4_str():
    """ Test four vector string representation. """
    vec = Vector4(1, 2, 3, 4)
    string = str(vec)
    assert string == 'FourVector(1.000000, 2.000000, 3.000000, 4.000000)'


def test_vec4_add():
    """ Test four vector addition. """
    vec1 = Vector4(1, 2, 3, 4)
    vec2 = Vector4(4, 3, 2, 1)
    vec3 = Vector4(5, 5, 5, 5)
    scalar = 3

    assert vec3 == vec1+vec2
    with pytest.raises(Exception):
        print(vec1+scalar)


def test_vec4_neg():
    """ Test four vector negation. """
    vec1 = Vector4(1, 2, 3, 4)
    vec2 = Vector4(-1, -2, -3, -4)

    assert vec1 == -vec2


def test_vec4_sub():
    """ Test four vector subtraction. """
    vec1 = Vector4(1, 2, 3, 4)
    vec2 = Vector4(1, 2, 3, 4)
    vec3 = Vector4()

    assert vec3 == vec1-vec2


def test_vec4_mul():
    """ Test four vector multiplication. """
    scalar1 = 3
    vec1 = Vector4(1, 2, 3, 4)
    vec2 = Vector4(3, 6, 9, 12)
    scalar2 = 16.0-9.0-4.0-1.0

    assert vec1*vec1 == scalar2
    assert vec1*vec1 == vec1.dot(vec1)
    assert scalar1*vec1 == vec1*scalar1
    assert scalar1*vec1 == vec2
    # assert vec1.__rmul__(vec2) == vec2*vec1
    with pytest.raises(Exception):
        print(vec1.dot(scalar2))


def test_vec4_div():
    """ Test four vector division. """
    scalar = 3.0
    vec1 = Vector4(1, 2, 3, 4)
    vec2 = Vector4(3, 6, 9, 12)

    assert vec2/scalar == vec1
    with pytest.raises(Exception):
        print(vec1/vec2)


def test_vec4_mass():
    """ Test four vector mass. """
    vec1 = Vector4(1., 2., 3., 4.)
    scalar = 2.

    assert vec1.mass2() == scalar
    assert vec1.mass() == np.sqrt(scalar)
    # assert vec1.mass2()**2 == vec1**4
    with pytest.raises(Exception):
        print(vec1**vec1)


def test_vec4_mom():
    """ Test four vector momentum functions. """
    vec1 = Vector4(1., 2., 3., 4.)
    pvec2 = 14.
    p_t2 = 5.

    assert vec1.p2() == pvec2
    assert vec1.p() == np.sqrt(pvec2)
    assert vec1.pt2() == p_t2
    assert vec1.pt() == np.sqrt(p_t2)


def test_vec4_angles():
    """ Test four vector angle functions. """
    vec = Vector4(0., 0., 4., 4.)
    vec1 = Vector4(1., 1., 1., 4.)
    vec2 = Vector4(1., -1., 1., 4.)

    assert vec.theta() == 0
    assert vec.phi() == 0
    assert vec1.phi() == np.pi/4.0
    assert vec2.phi() == 7.0*np.pi/4.0


def test_vec4_cross():
    """ Test four vector cross product. """
    vec1 = Vector4(3., 2., 1., 4.0)
    vec2 = Vector4(1., 2., 3., 4.0)
    vec3 = Vector4(4., -8., 4., 0.0)
    scalar = 5.

    assert vec1.cross(vec2) == vec3
    assert vec1.cross(vec2) == -vec2.cross(vec1)

    with pytest.raises(Exception):
        print(vec1.cross(scalar))


def test_vec4_boost():
    """ Test four vector boost. """
    vec1 = Vector4(3., 2., 1., 4.)
    vec2 = Vector4(12., 2., -3., 25.)
    beta = vec2.boost_vector()
    vec4 = vec1.boost(beta).boost(-beta)

    assert abs(vec1.px()-vec4.px()) < 1E-8
    assert abs(vec1.py()-vec4.py()) < 1E-8
    assert abs(vec1.pz()-vec4.pz()) < 1E-8
    assert abs(vec1.energy()-vec4.energy()) < 1E-8
    with pytest.raises(Exception):
        print(vec1.boost(vec2))
