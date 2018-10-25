from .context import nuChic
from nuChic.FourVector import Vec4
from nuChic.ThreeVector import Vec3
import numpy as np
import pytest

def test_Vec4_init():
    v1 = Vec4(1,2,3,4)
    v2 = Vec4()
    v2.E = 1
    v2.px = 2
    v2.py = 3
    v2.pz = 4
    assert v1 == v2

def test_Vec4_getitem():
    v = Vec4(1,2,3,4)
    assert v[0] == 1
    assert v[1] == 2
    assert v[2] == 3
    assert v[3] == 4
    with pytest.raises(Exception):
        v[4]

def test_Vec4_repr():
    v = Vec4(1,2,3,4)
    string = repr(v)
    v2 = eval(string)
    assert v == v2

def test_Vec4_str():
    v = Vec4(1,2,3,4)
    string = str(v)
    assert string == '(1,2,3,4)'

def test_Vec4_add():
    v1 = Vec4(1,2,3,4)
    v2 = Vec4(4,3,2,1)
    v3 = Vec4(5,5,5,5)
    b = 3

    assert v3 == v1+v2
    with pytest.raises(Exception):
        v1+b

def test_Vec4_neg():
    v1 = Vec4(1,2,3,4)
    v2 = Vec4(-1,-2,-3,-4)

    assert v1 == -v2

def test_Vec4_sub():
    v1 = Vec4(1,2,3,4)
    v2 = Vec4(1,2,3,4)
    v3 = Vec4()

    assert v3 == v1-v2

def test_Vec4_mul():
    a = 3
    v1 = Vec4(1,2,3,4)
    v2 = Vec4(3,6,9,12)
    b = -28

    assert v1*v1 == b
    assert v1*v1 == v1.dot(v1)
    assert a*v1 == v1*a
    assert a*v1 == v2
    with pytest.raises(Exception):
        v1.dot(b)

def test_Vec4_div():
    a = 3.0
    v1 = Vec4(1,2,3,4)
    v2 = Vec4(3,6,9,12)

    assert v2/a == v1
    with pytest.raises(Exception):
        v1/v2

def test_Vec4_Mass():
    v1 = Vec4(4,3,2,1)
    b = 2

    assert v1.M2() == b
    assert v1.M() == np.sqrt(b)

def test_Vec4_P():
    v1 = Vec4(4,3,2,1)
    pvec2 = 14
    pt2 = 13

    assert v1.P2() == pvec2
    assert v1.P() == np.sqrt(pvec2)
    assert v1.Pt2() == pt2
    assert v1.Pt() == np.sqrt(pt2)

def test_Vec4_Angles():
    v = Vec4(4,0,0,4)
    v1 = Vec4(4,1,1,1)

    assert v.Theta() == 0
    assert v.Phi() == 0
    assert v1.Phi() == np.pi/4.0

def test_Vec4_Cross():
    v1 = Vec4(4,3,2,1)
    v2 = Vec4(4,1,2,3)
    v3 = Vec4(0.0,4,-8,4)
    a = 5 

    assert v1.Cross(v2) == v3
    assert v1.Cross(v2) == -v2.Cross(v1) 

    with pytest.raises(Exception):
        v1.Cross(a)

def test_Vec4_Boost():
    v1 = Vec4(4,3,2,1)
    v2 = Vec4(25,12,2,-3)
    beta = v2.BoostVector()
    v4 = v1.Boost(beta).BoostBack(beta)

    assert abs(v1.px-v4.px) < 1E-8
    assert abs(v1.py-v4.py) < 1E-8
    assert abs(v1.pz-v4.pz) < 1E-8
    assert abs(v1.E-v4.E) < 1E-8
