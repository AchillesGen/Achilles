import pandas as pd
from scipy import interpolate
from scipy.signal import savgol_filter
import numpy as np

from nuChic.Particle import Particle
from nuChic.Nucleus import Nucleus
from nuChic.Cascade import FSI
from nuChic.Constants import mN, MeV, GeV
from nuChic.FourVector import Vec4
from nuChic.ThreeVector import Vec3

import matplotlib.pyplot as plt

c12Density_db = pd.read_csv('densities/c12.density',header=10,sep='\s+',names=['r','rho','error'])
c12Density = interpolate.InterpolatedUnivariateSpline(c12Density_db['r'].values,c12Density_db['rho'],k=1)

r = np.linspace(0.,5,1000)
rho = c12Density(r)

density_smooth = savgol_filter(c12Density_db['r'].values, 3, 1) # window size 3, polynomial order 1
c12Density_smooth = interpolate.InterpolatedUnivariateSpline(density_smooth,c12Density_db['rho'])
rho_smooth = c12Density_smooth(r)

C12 = Nucleus(6,12,92.15,225,c12Density,c12Density)
energy = np.linspace(C12.kf * MeV, 3000 * MeV, 100)
mean_escape = np.zeros_like(energy)

for j, p in enumerate(energy):
    phi = 2.0*np.pi*np.random.uniform()
    costheta = 2.0*np.random.uniform() - 1.0
    sintheta = np.sqrt(1.0 - costheta**2)
    px = p * sintheta * np.sin(phi)
    py = p * sintheta * np.cos(phi)
    pz = p * costheta
    
    prop_dist = []
    nevents = 1000
    for i in range(nevents):
        particle = Particle(2112, Vec4(np.sqrt(mN**2 + px**2 + py**2 + pz**2), px, py, pz), Vec3(0,0,0))
        fsi = FSI(C12, 3e-1)
        dist = fsi([particle])
        if dist is not None:
            prop_dist.append(len(dist))

    mean_escape[j] = np.mean(prop_dist)
    print(p, mean_escape[j])

plt.plot(energy, mean_escape)
plt.show()


