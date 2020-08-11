import numpy as np
import matplotlib.pyplot as plt
import nuchic.constants as constants
import nuchic.physics as physics
import nuchic.utilities as utilities


class Oset:
    def __init__(self, Z, N, fpi = np.sqrt(0.360*4*np.pi)):
        self.chi = (N-Z)/(N+Z)
        self.fpi = fpi
        self.pmat = np.array([[5-4*self.chi, 1-self.chi, 0],
                              [1+self.chi, 4, 1-self.chi],
                              [0, 1+self.chi, 5+4*self.chi]])
        self.m_delta = physics.ParticleInfo(2114).mass()
        self.m_pi = physics.ParticleInfo(111).mass()
        self.m_n = physics.ParticleInfo(2212).mass()

    def __call__(self, radius, ppi):
        swave = self.swave(radius, ppi)
        pwave = self.pwave(radius, ppi)
        abs_swave = self.absorbtion_swave(ppi.energy(), radius)
        abs_pwave = self.absorbtion_pwave(radius, ppi)
        return (pwave+swave)/197.3269804, (abs_pwave+abs_swave)/197.3269804

    def pwave(self, radius, ppi):
        den = self.density(radius)*(197.3269804)**3
        kf = (3/2*np.pi**2*den)**(1/3)
        pn = physics.Vector4(0, 0, kf, np.sqrt(self.m_n**2+kf**2))
        ptot = ppi+pn
        ecm = np.sqrt(2*ppi.energy()*self.m_n)
        qcm = ppi.pz()*self.m_n/ecm
        gamma = 1.0/(12.0*np.pi)*(self.fpi/m_pi)**2*self.m_n/ecm*qcm**3
        prop = 1/(ecm - self.m_delta + 1j*gamma)
        prop = (prop*prop.conjugate()).real
        probability = 1/ppi.energy()*den*2/3*(self.fpi/self.m_pi)**2*qcm**2*prop*gamma
        return probability/3*(2-self.chi)


    def swave(self, radius, ppi):
        den = self.density(radius)*(197.3269804)**3
        kf = (3/2*np.pi**2*den)**(1/3)
        pn = physics.Vector4(0, 0, kf, np.sqrt(self.m_n**2+kf**2))
        ptot = ppi+pn
        ecm = ptot.mass()
        xi = (ecm - self.m_n - self.m_pi)/self.m_pi
        sigma = (0.19753+0.06899*xi-0.01334*xi**2)/self.m_pi**2
        D = -0.03130+0.37062*xi-0.08229*xi**2
        B = 0.21972+0.06602*xi-0.01866*xi**2 
        C = 0.5*(1-D)
        return ppi.pz()/ppi.energy()*sigma*den*(1-(B-C)*self.chi)
    
    def density(self, radius):
        return 0.16/(1+np.exp((radius-3.971)/0.5935))

    def absorbtion_swave(self, energy, radius):
        den = self.density(radius)*(197.3269804)**3
        return 4*np.pi/energy*(1+energy/(2*self.m_n))*0.035*self.m_pi**(-4)*den**2

    def absorbtion_pwave(self, radius, ppi):
        den = self.density(radius)*(197.3269804)**3
        den_0 = self.density(0)*(197.3269804)**3
        CQ, CA2, CA3, alpha, beta, gamma = 12.0, 16.3, 15.8, 0.42, 0.80, 1.60
        sigma_delta = -(CQ*(den/den_0)**alpha + CA2*(den/den_0)**beta + CA3*(den/den_0)**gamma)
        sigma_delta2 = -(CA2*(den/den_0)**beta + CA3*(den/den_0)**gamma)
        kf = (3/2*np.pi**2*den)**(1/3)
        pn = physics.Vector4(0, 0, kf, np.sqrt(self.m_n**2+kf**2))
        ptot = ppi+pn
        ecm = np.sqrt(ppi.energy()*self.m_n)
        qcm = ppi.pz()*self.m_n/ecm
        gamma = 1.0/(12.0*np.pi)*(self.fpi/m_pi)**2*self.m_n/ecm*qcm**3
        prop = 1/(ecm - self.m_delta + 1j*(gamma-sigma_delta))
        prop = (prop*prop.conjugate()).real
        return 1/ppi.energy()*4/9*(self.fpi/self.m_pi)**2*qcm**2*prop*(gamma-sigma_delta2)*den



if __name__ == '__main__':
    m_pi = physics.ParticleInfo(111).mass()
    m_n = physics.ParticleInfo(2212).mass()
    tpi = 165
    epi = tpi+m_pi
    ppi = physics.Vector4(0, 0, np.sqrt(epi**2-m_pi**2), epi)

    oset = Oset(26, 30)

    radii = np.linspace(0, 7, 100)
    rho = oset.density(radii)*(197.3269804/m_pi)**3
    scatter = []
    absorbed = []
    for radius in radii:
        interact, absorb = oset(radius, ppi)
        scatter.append(interact)
        absorbed.append(absorb)
    print(scatter)
    print(absorbed)
    plt.plot(radii, rho)
    plt.plot(radii, scatter)
    plt.plot(radii, absorbed)
    plt.xlim([0,7])
    plt.ylim([0,0.8])
    plt.show()
