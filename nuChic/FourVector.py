import numpy as np

class Vec4:
    def __init__(self, E=0, px=0, py=0, pz=0):
        self.E = E
        self.px = px
        self.py = py
        self.pz = pz

    def __getitem__(self,i):
        if i == 0: return self.E
        if i == 1: return self.px
        if i == 2: return self.py
        if i == 3: return self.pz
        raise Exception('Vec4D')

    def __repr__(self):
        return 'Vec4({0},{1},{2},{3})'.format(self.E, self.px, self.py, self.pz)

    def __str__(self):
        return '({0},{1},{2},{3})'.format(self.E, self.px, self.py, self.pz)
    
    def __add__(self,v):
        if not isinstance(v,Vec4):
            raise Exception('Vec4')
        return Vec4(self.E+v.E, self.px+v.px, self.py+v.py, self.pz+v.pz)

    def __neg__(self):
        return Vec4(-self.E,-self.px,-self.py,-self.pz)

    def __sub__(self,v):
        return self + -v

    def __mul__(self,v):
        if isinstance(v,Vec4):
            return self.dot(v)
        return Vec4(self.E*v,self.px*v,self.py*v,self.pz*v)

    def __rmul__(self,v):
        if isinstance(v,Vec4):
            return self.dot(v)
        return Vec4(self.E*v,self.px*v,self.py*v,self.pz*v)

    def __truediv__(self,v):
        if isinstance(v,Vec4):
            raise Exception('Vec4')
        return (1.0/v)*self

    def __eq__(self,v):
        if isinstance(v,Vec4):
            if(self.E == v.E and
               self.px == v.px and
               self.py == v.py and
               self.pz == v.pz):
                return True
        return False

    def dot(self,v):
        if not isinstance(v,Vec4):
            raise Exception('Vec4')
        return self.E*v.E - (self.px*v.px + self.py*v.py + self.pz*v.pz)

    def M2(self):
        return self.dot(self)

    def M(self):
        return np.sqrt(self.M2())

    def P2(self):
        return self.px*self.px+self.py*self.py+self.pz*self.pz

    def P(self):
        return np.sqrt(self.P2())

    def Pt2(self):
        return self.px*self.px+self.py*self.py

    def Pt(self):
        return np.sqrt(self.Pt2())

    def Theta(self):
        return np.arccos(self.pz/self.P())

    def Phi(self):
        if self.px == 0 and self.py == 0:
            return 0.0
        else:
            return np.arctan2(self.py,self.px)

    def Cross(self,v):
        if not isinstance(v,Vec4):
            raise Exception('Vec4')
        return Vec4(0.0,
                    self.py*v.pz - self.pz*v.py,
                    self.pz*v.px - self.px*v.pz,
                    self.px*v.py - self.py*v.px)

    def Boost(self,v):
        if not isinstance(v,Vec4):
            raise Exception('Vec4')

        beta = (v.px/v.E,v.py/v.E,v.pz/v.E)
        beta2 = sum(x*x for x in beta)
        gamma = 1.0/np.sqrt(1.0-beta2)
        betap = beta[0]*self.px+beta[1]*self.py+beta[2]*self.pz
        gamma2 = (gamma-1.0)/beta2 if beta2 > 0 else 0.0

        return Vec4(gamma*(self.E+betap),
            self.px+gamma2*betap*beta[0]+gamma*beta[0]*self.E,
            self.py+gamma2*betap*beta[1]+gamma*beta[1]*self.E,
            self.pz+gamma2*betap*beta[2]+gamma*beta[2]*self.E)

    def BoostBack(self,v):
        if not isinstance(v,Vec4):
            raise Exception('Vec4')

        beta = (-v.px/v.E,-v.py/v.E,-v.pz/v.E)
        beta2 = sum(x*x for x in beta)
        gamma = 1.0/np.sqrt(1.0-beta2)
        betap = beta[0]*self.px+beta[1]*self.py+beta[2]*self.pz
        gamma2 = (gamma-1.0)/beta2 if beta2 > 0 else 0.0

        return Vec4(gamma*(self.E+betap),
            self.px+gamma2*betap*beta[0]+gamma*beta[0]*self.E,
            self.py+gamma2*betap*beta[1]+gamma*beta[1]*self.E,
            self.pz+gamma2*betap*beta[2]+gamma*beta[2]*self.E)

