#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from nuChic.Nucleus import Nucleus
from nuChic.Constants import hbarc, MeV, GeV, fm, mN
from nuChic.Cascade import FSI
import time
import pandas as pd


argon_nucleus = Nucleus(6,12, 8.6*MeV, 225*MeV)
start_time = time.time()
points = 1000
number_of_outgoing_particles =[]
foo = FSI(argon_nucleus, 500.*MeV, 1)
for i in range(points):
    foo.kick(500.*MeV)
    number_of_outgoing_particles.append(len(foo()))
#    number_of_outgoing_particles.append(len(foo.outgoing_particles))
#    print("Number of outgoing nucleons: ",len(foo.outgoing_particles))
    foo.reset()
print("Done!")
print("--- %s seconds ---" % (time.time() - start_time))
print("--- %s seconds/point ---" % ((time.time() - start_time)/points))
print("Average number of outgoing nucleons: ",sum(number_of_outgoing_particles)/points)

#print("Average number of outgoing nucleons: ")
##line = linecache.getline("/Users/pmachado/Dropbox/Projects/NuGen/FNALNeuGen/configurations/pos_part_in_v2.out", 33)
#c12Density_db = pd.read_csv("/Users/pmachado/Dropbox/Projects/NuGen/FNALNeuGen/configurations/pos_part_in_v2.out",sep='\s+',names=['index','pid','x','y','z'], nrows=100)
#
##print(c12Density_db)
##temp = np.asarray(c12Density_db.iloc[2:5])
#i0=12
#protons = np.asarray(c12Density_db.iloc[i0:i0+6][['x','y','z']])
#neutrons = np.asarray(c12Density_db.iloc[i0+6:i0+12][['x','y','z']])
#
#print('PROTONS')
#print(protons)
#print('NEUTRONS')
#print(neutrons)
#print('done')
