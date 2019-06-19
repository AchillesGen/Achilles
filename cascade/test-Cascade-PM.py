#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from nuChic.Nucleus import Nucleus
from nuChic.Constants import hbarc, MeV, GeV, fm, mN
from CascadePM import FSI


argon_nucleus = Nucleus(18,40, 8.6*MeV, 200*MeV)
for i in range(50):
    foo = FSI(argon_nucleus, 100.*MeV, 1)
    foo()
print("Done!")
print(foo.outgoing_particles[0].mom)
print("Number of outgoing nucleons: ",len(foo.outgoing_particles))