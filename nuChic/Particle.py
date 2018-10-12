import numpy as np
from FourVector import Vec4

class Particle:
    def __init__(self,pid=0,mom=Vec4(),status=0,mothers=[],daughters=[]):
        self.pid = pid
        self.mom = mom
        self.status = status
        self.mothers = mothers
        self.daughters = daughters

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
