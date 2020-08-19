import numpy as np
import matplotlib.pyplot as plt
import nuchic.constants as constants
import nuchic.physics as physics
import nuchic.utilities as utilities


HBARC = constants.HBARC


class Oset:
    def __init__(self, Z, N, fstr2 = 0.360*4*np.pi):
        self.chi = (N-Z)/(N+Z)
        self.pmat = np.array([[5-4*self.chi, 1-self.chi, 0],
                              [1+self.chi, 4, 1-self.chi],
                              [0, 1+self.chi, 5+4*self.chi]])
        self.m_delta = physics.ParticleInfo(2114).mass()
        self.m_pi = (physics.ParticleInfo(211).mass()+physics.ParticleInfo(111).mass())/2.0
        self.m_n = physics.ParticleInfo(2212).mass()
        self.norm = HBARC**2 # conversion from MeV^-2 to fm^2
        self.coupling = fstr2/self.m_pi**2
        self.rho0 = 0.16843764

    def pwave(self, costh, mom):
        sinth = np.sqrt(1-costh**2)
        k = physics.Vector4(mom*sinth, 0, mom*costh, np.sqrt(self.m_n**2+mom**2))

        pcm4 = ppi.boost(-(k+ppi).boost_vector())

        # ecm = (k+self.ppi).mass()
        # pcm = pcm4.p()
        ecm = np.sqrt(self.m_pi**2 + 2*self.m_n*ppi.energy() + self.m_n**2)
        pcm = ppi.p()*self.m_n/ecm
        self.self_energy_calc(self.rho)

        gamma = 1.0/(12.0*np.pi)*self.coupling*self.m_n/ecm*pcm**3*self.pauli_blocking(self.ppi, ecm, pcm, self.kf, k)
        prop2 = 1/((ecm - self.m_delta)**2 + (gamma-self.self_energy)**2)

        pqel = 1/ppi.energy()*2*mom**2/(np.pi**2)/3*self.coupling*pcm**2*prop2*gamma/HBARC

        return  pqel


    def __call__(self, radius, ppi):
        self.rho = self.density(radius)
        self.kf = (3/2*np.pi**2*self.rho)**(1/3)*HBARC
        self.ppi = ppi
        result = dblquad(self.pwave, 0, self.kf, -1, 1)

        ecm = np.sqrt(self.m_pi**2 + 2*self.m_n*ppi.energy() + self.m_n**2)
        xi = (ecm - self.m_n - self.m_pi)/self.m_pi
        sqel = (0.19753 + 0.06899*xi - 0.01334*xi**2)/self.m_pi**2*self.norm*self.rho*self.ppi.p()/self.ppi.energy()
        print(result, sqel)
        return result[0]+sqel
        
        # self.pxsec = self.coupling*prop2*pcm**2*2/3
        # 
        # # Absorption
        # abs_pwave = 4.0/9.0*self.pxsec * (-self.self_energy_abs)
        # imb0 = 0.035/(self.m_pi)**4
        # abs_swave = 4*np.pi/ppi.energy()*rho**2*(1+ppi.energy()/(2*self.m_n))*imb0*HBARC**5

        # # Qel
        # pqel = self.pxsec/ppi.energy() * 4 * gamma * mom * mom / (2*np.pi)**3 * 4*np.pi * kf / HBARC
        # print(pqel, mom * mom * (4*np.pi)/(2*np.pi)**3 * kf, rho)
        # D = -0.03130 + 0.37062*xi - 0.08229*xi**2
        # B = 0.21972 + 0.06602*xi - 0.01866*xi**2
        # A = 0.5*(1+D)
        # C = 0.5*(1-D)
        # pqel_tot = pqel*(2-self.chi)/3
        # sqel_tot = sqel*(1-(B-C)*self.chi)

        # prob_factor = rho*ppi.p()/ppi.energy()
        # return (pqel_tot+sqel_tot)*prob_factor, abs_pwave*prob_factor+abs_swave

    def self_energy_calc(self, rho):
        rho_ratio = rho/self.rho0
        alpha = 0.42 
        beta = 0.80
        gamma1 = 1.60
        CQ = 12.0
        CA2 = 16.3
        CA3 = 15.8
        abs2 = CA2*(rho_ratio)**beta
        abs3 = CA3*(rho_ratio)**gamma1
        self.self_energy_abs = - abs2 - abs3
        self.self_energy = self.self_energy_abs - CQ*(rho_ratio)**alpha

    def pauli_blocking(self, ppi, ecm, pcm, kf, k):
        ef = np.sqrt(kf**2+self.m_n**2)
        edelta = ppi.energy()+self.m_n
        edcm = np.sqrt(pcm*pcm + self.m_n**2)
        mu0 = (edelta*edcm - ef*ecm)/(ppi+k).p()/pcm
        if(mu0 < -1):
            return -1
        elif(mu0 > 1):
            return 1
        return (mu0**3+mu0+2)/4

    def density(self, radius):
        return self.rho0/(1+np.exp((radius-3.93341372)/0.55863863))


def fit(x, a, b, c):
    return c/(1+np.exp((x-a)/b))*(HBARC/m_pi)**3


def test(k):
    return 4 * k*k/(2*np.pi)**3 * (4*np.pi)


if __name__ == '__main__':
    import pandas as pd
    from scipy.optimize import curve_fit
    from scipy.integrate import quad, dblquad

    result = quad(test, 0, 225)
    print((3*np.pi**2/2*result[0])**(1/3))

    m_pi = (physics.ParticleInfo(211).mass()+physics.ParticleInfo(111).mass())/2.0
    m_n = physics.ParticleInfo(2212).mass()
    tpi = 165
    epi = tpi+m_pi
    ppi = physics.Vector4(0, 0, np.sqrt(epi**2-m_pi**2), epi)

    oset = Oset(26, 30)

    radii = np.linspace(0, 7, 50)
    rho = oset.density(radii)*(HBARC/m_pi)**3
    scatter = []
    absorbed = []
    for radius in radii:
        interact = oset(radius, ppi)
        scatter.append(interact)

    oset_data = pd.read_csv('oset.csv',
                            names=['Qx', 'Qy', 'Ax', 'Ay', 'rhox', 'rhoy'])

    popt, pcov = curve_fit(fit, oset_data['rhox'][:46], oset_data['rhoy'][:46])

    # Plot Oset data
    plt.plot(radii, scatter, color='tab:red')
    # plt.plot(radii, absorbed, color='tab:blue')
    plt.plot(radii, rho, color='tab:green')

    # Plot our calculation
    plt.plot(oset_data['Qx'], oset_data['Qy'], color='tab:red', ls='--')
    plt.plot(oset_data['Ax'], oset_data['Ay'], color='tab:blue', ls='--')
    plt.plot(oset_data['rhox'], oset_data['rhoy'], color='tab:green', ls='--')

    plt.xlim([0,7])
    plt.ylim([0,0.8])
    plt.show()
