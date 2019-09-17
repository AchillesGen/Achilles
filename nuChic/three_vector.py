"""Spatial three-vector"""
import numpy as np


class Vec3:
    """Spatial three-vector class.

    Attributes
        x: the x-component
        y: the y-component
        z: the z-component
    """

    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __getitem__(self, i):
        if i == 0:
            return self.x
        if i == 1:
            return self.y
        if i == 2:
            return self.z
        raise Exception('Vec3D')

    def __repr__(self):
        return 'Vec3({0},{1},{2})'.format(self.x, self.y, self.z)

    def __str__(self):
        return '({0},{1},{2})'.format(self.x, self.y, self.z)

    def __add__(self, v):
        if not isinstance(v, Vec3):
            raise Exception('Vec3D')
        return Vec3(self.x+v.x, self.y+v.y, self.z+v.z)

    def __neg__(self):
        return Vec3(-self.x, -self.y, -self.z)

    def __sub__(self, v):
        return self + -v

    def __mul__(self, v):
        if isinstance(v, Vec3):
            return self.dot(v)
        return Vec3(self.x*v, self.y*v, self.z*v)

    def __rmul__(self, v):
        if isinstance(v, Vec3):
            return self.dot(v)
        return Vec3(self.x*v, self.y*v, self.z*v)

    def __truediv__(self, v):
        if isinstance(v, Vec3):
            raise Exception('Vec3')
        return (1.0/v)*self

    def __eq__(self, v):
        if isinstance(v, Vec3):
            if(self.x == v.x and
               self.y == v.y and
               self.z == v.z):
                return True
        return False

    @property
    def vec(self):
        """ The three-vector as a list """
        return [self.x, self.y, self.z]

    @property
    def array(self):
        """ The three-vector as a np.array """
        return np.array([self.x, self.y, self.z])

    def dot(self, v):
        """ The dot product of the vector with v """
        if not isinstance(v, Vec3):
            raise Exception('Vec3')
        return self.x*v.x + self.y*v.y + self.z*v.z

    @property
    def p2(self):
        """ The square of the vector """
        return self.dot(self)

    @property
    def p(self):
        """ The modulus of the vector """
        return np.sqrt(self.p2)

    @property
    def theta(self):
        """ The polar angle theta of spherical coordinates """
        return np.arccos(self.z/self.p)

    @property
    def phi(self):
        """ The azimuthal angle theta of spherical coordinates """
        if self.x == 0 and self.y == 0:
            return 0.0
        return np.arctan2(self.y, self.x)

    def cross(self, vec):
        """ The cross product of three-vectors (p x v)_i = eps_{ijk} p_j v_k"""
        if not isinstance(vec, Vec3):
            raise Exception('Vec3')
        return Vec3(self.y*vec.z - self.z*vec.y,
                    self.z*vec.x - self.x*vec.z,
                    self.x*vec.y - self.y*vec.x)
