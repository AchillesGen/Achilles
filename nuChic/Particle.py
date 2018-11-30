import numpy as np
from nuChic.FourVector import Vec4
from nuChic.ThreeVector import Vec3

class Particle:
    def __init__(self,pid=0,mom=Vec4(),pos=Vec3(),status=0,mothers=None,daughters=None):
        self.pid = pid
        self.mom = mom
        self.pos = pos
        self.status = status
        self.mothers = mothers if mothers is not None else []
        self.daughters = daughters if daughters is not None else []
        self.vec = Vec3(self.mom[1]/self.mom[0],self.mom[2]/self.mom[0],self.mom[3]/self.mom[0])

    def isFinal(self):
        if self.status == 1:
            return True
        return False

    def Px(self):
        return self.mom.px

    def Py(self):
        return self.mom.py

    def Pz(self):
        return self.mom.pz

    def E(self):
        return self.mom.E

    def M(self):
        return self.mom.M()

    def r(self):
        return self.pos.P()

    def propagate(self,path,time):
        # Propagation of a particle with mean free path=path with a time step(time)
        dist = self.vec.P()*time

        if path < dist:
            prop_dist = path
        else:
            prop_dist = dist

        theta = self.vec.Theta()
        phi = self.vec.Phi()

        prop_dist_x = prop_dist*np.sin(theta)*np.cos(phi)
        prop_dist_y = prop_dist*np.sin(theta)*np.sin(phi)
        prop_dist_z = prop_dist*np.cos(theta)

        self.prop_dist = Vec3(prop_dist_x,prop_dist_y,prop_dist_z)
        self.pos += self.prop_dist

        return prop_dist

    def back_propagate(self):
        self.pos -= self.prop_dist
