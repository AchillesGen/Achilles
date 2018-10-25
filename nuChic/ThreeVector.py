import numpy as np

class Vec3:
    def __init__(self,x=0,y=0,z=0):
        self.x = x
        self.y = y
        self.z = z

    def __getitem__(self,i):
        if i == 0: return self.x
        if i == 1: return self.y
        if i == 2: return self.z
        raise Exception('Vec3D')

    def __repr__(self):
        return 'Vec3({0},{1},{2})'.format(self.x,self.y,self.z)

    def __str__(self):
        return '({0},{1},{2})'.format(self.x,self.y,self.z)

    def __add__(self,v):
        if not isinstance(v,Vec3):
            raise Exception('Vec3D')
        return Vec3(self.x+v.x, self.y+v.y, self.z+v.z)

    def __neg__(self):
        return Vec3(-self.x,-self.y,-self.z)

    def __sub__(self,v):
        return self + -v

    def __mul__(self,v):
        if isinstance(v,Vec3):
            return self.dot(v)
        return Vec3(self.x*v, self.y*v, self.z*v)

    def __rmul__(self,v):
        if isinstance(v,Vec3):
            return self.dot(v)
        return Vec3(self.x*v, self.y*v, self.z*v)

    def __truediv__(self,v):
        if isinstance(v,Vec3):
            raise Exception('Vec3')
        return (1.0/v)*self

    def __eq__(self,v):
        if isinstance(v,Vec3):
            if(self.x == v.x and
               self.y == v.y and
               self.z == v.z):
                return True
        return False

    def Vec(self):
        return [self.x, self.y, self.z]

    def dot(self,v):
        if not isinstance(v,Vec3):
            raise Exception('Vec3')
        return self.x*v.x + self.y*v.y + self.z*v.z

    def P2(self):
        return self.dot(self)

    def P(self):
        return np.sqrt(self.P2())

    def Theta(self):
        return np.arccos(self.z/self.P()) 

    def Phi(self):
        if self.x == 0 and self.y == 0:
            return 0.0
        else:
            return np.arctan2(self.y,self.x)

    def Cross(self,v):
        if not isinstance(v,Vec3):
            raise Exception('Vec3')
        return Vec3(self.y*v.z - self.z*v.y,
                    self.z*v.x - self.x*v.z,
                    self.x*v.y - self.y*v.x)
