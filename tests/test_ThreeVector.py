from .context import nuChic
from nuChic.ThreeVector import Vec3
import numpy as np
import pytest

def test_Vec3_init():
    v1 = Vec3(1,2,3)
    v2 = Vec3()
    v2.x = 1
    v2.y = 2
    v2.z = 3
    v3 = Vec3(3,2,1)

    assert v1 == v2
    assert v1 != v3

def test_Vec3_getitem():
    v = Vec3(1,2,3)
    assert v[0] == 1
    assert v[1] == 2
    assert v[2] == 3
    with pytest.raises(Exception):
        v[3]

def test_Vec3_repr():
    v = Vec3(1,2,3)
    string = repr(v)
    v2 = eval(string)
    assert v == v2

def test_Vec3_str():
    v = Vec3(1,2,3)
    string = str(v)
    assert string == '(1,2,3)'

def test_Vec3_add():
    v1 = Vec3(1,2,3)
    v2 = Vec3(4,3,2)
    v3 = Vec3(5,5,5)
    b = 3

    assert v3 == v1+v2
    with pytest.raises(Exception):
        v1+b

def test_Vec3_neg():
    v1 = Vec3(1,2,3)
    v2 = Vec3(-1,-2,-3)

    assert v1 == -v2

def test_Vec3_sub():
    v1 = Vec3(1,2,3)
    v2 = Vec3(1,2,3)
    v3 = Vec3()

    assert v3 == v1-v2

def test_Vec3_mul():
    a = 3
    v1 = Vec3(1,2,3)
    v2 = Vec3(3,6,9)
    b = 14

    assert v1*v1 == b
    assert v1*v1 == v1.dot(v1)
    assert a*v1 == v1*a
    assert a*v1 == v2
    assert v1.__rmul__(v2) == v2*v1

    with pytest.raises(Exception):
        v1.dot(b)

def test_Vec3_div():
    a = 3.0
    v1 = Vec3(1,2,3)
    v2 = Vec3(3,6,9)

    assert v2/a == v1
    with pytest.raises(Exception):
        v1/v2

def test_Vec3_P():
    v1 = Vec3(4,3,2)
    pvec2 = 29

    assert v1.P2() == pvec2
    assert v1.P() == np.sqrt(pvec2)

def test_Vec3_Angles():
    v = Vec3(0,0,4)
    v1 = Vec3(1,1,1)

    assert v.Theta() == 0
    assert v.Phi() == 0
    assert v1.Phi() == np.pi/4.0

def test_Vec3_Cross():
    v1 = Vec3(3,2,1)
    v2 = Vec3(1,2,3)
    v3 = Vec3(4,-8,4)
    a = 5 

    assert v1.Cross(v2) == v3
    assert v1.Cross(v2) == -v2.Cross(v1) 

    with pytest.raises(Exception):
        v1.Cross(a)
