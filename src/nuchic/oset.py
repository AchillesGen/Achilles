import numpy as np
import matplotlib.pyplot as plt
import nuchic.constants as constants
import nuchic.physics as physics
import nuchic.utilities as utilities


HBARC = constants.HBARC


class Oset:
    def __init__(self, Z, N, fstr2 = 0.360*4*np.pi):
        self.chi = (N-Z)/(N+Z)
        self.Z = Z
        self.N = N
        self.A = Z+N
        self.pmat = np.array([[5-4*self.chi, 1-self.chi, 0],
                              [1+self.chi, 4, 1-self.chi],
                              [0, 1+self.chi, 5+4*self.chi]])
        self.m_delta = physics.ParticleInfo(2114).mass()
        self.m_pi = physics.ParticleInfo(211).mass()
        self.m_n = physics.ParticleInfo(2212).mass()
        self.norm = HBARC**2 # conversion from MeV^-2 to fm^2
        self.coupling = fstr2/self.m_pi**2
        self.rho0 = 0.17
        self.imb0 = 0.035/self.m_pi**4
        self.ax3 = 1/(12*np.pi)*self.coupling*self.m_n
        self.ax2 = 4/9*self.coupling

    def cq(self, x):
        return -5.19*x**2 + 15.35*x + 2.06

    def ca2(self, x):
        return 1.06*x**2 - 6.64*x + 22.66

    def ca3(self, x):
        return -13.46*x**2 + 46.17*x - 20.34

    def alpha(self, x):
        return 0.382*x**2 - 1.322*x + 1.466

    def beta(self, x):
        return -0.038*x**2 + 0.204*x + 0.613

    def __call__(self, cos_piN,kf, radius, ppi,flag):
        rhop=self.density(radius)
        rhon=self.density(radius)/self.Z*self.N

        rho = (rhop+rhon)*0.5
        #kfp = (3*np.pi**2*rhop)**(1/3)*HBARC
        #kfn = (3*np.pi**2*rhon)**(1/3)*HBARC
        #kf=0.5*(kfp+kfn)
        ef = np.sqrt(kf**2 + self.m_n**2)

        # Absorption
        wqa = ppi.energy()
        if(ppi.energy() - ppi.mass() > 0.5*HBARC):
            wqa = self.m_pi + 0.5*HBARC 
        # Eq. 3.13         
        #abs_swave = 4*np.pi/ppi.p()*(1+wqa/(2*self.m_n))*self.imb0*rho*rho
        abs_swave = 4*np.pi/wqa*(1+wqa/(2*self.m_n))*self.imb0*rhop*rhop
        # Convert to fm^-1
        abs_swave *= HBARC**5

        # Delta kinematics
        #p2 = ppi.p2() + 0.6*kf**2 +2.0*ppi.p()*np.sqrt(0.6)*kf*cos_piN
        p2 = ppi.p2() + kf**2 +2.0*ppi.p()*kf*cos_piN
        p = np.sqrt(p2)
        edelta = ppi.energy() + np.sqrt(kf**2+self.m_n**2)
        s = edelta*edelta - p2
        sqrts = np.sqrt(s)

        # Pion CM
        wcm = (s-self.m_n**2+self.m_pi**2)/(2*sqrts)
        pcm2 = wcm*wcm - self.m_pi**2
        pcm = np.sqrt(pcm2)

        # Nucleon CM energy
        en = sqrts - wcm
        # eq 2.13
        mu0 = (edelta*en - ef*sqrts)/(p*pcm)
        if(mu0 < -1):
            mu0 = -1
        elif(mu0 > 1):
            mu0 = 1

        pauli = 0.25*(2+mu0+mu0**3)
        if(pauli < 0):
            pauli = 0
        elif(pauli > 1):
            pauli = 1

        gamma_free = self.ax3*pcm*pcm2/sqrts
        gamma = gamma_free*pauli
        x = (ppi.energy()-ppi.mass())/self.m_pi
        # eq 2.21
        rho_ratio = rho/self.rho0
        alpha = self.alpha(x)
        gamma_cq = self.cq(x)*rho_ratio**alpha
        if(gamma_cq < 0):
            gamma_cq = 0
        beta = self.beta(x)
        gamma_ca2 = self.ca2(x)*rho_ratio**beta
        if(gamma_ca2 < 0):
            gamma_ca2 = 0
        gamma_ca3 = self.ca3(x)*rho_ratio**(2*beta)
        if(gamma_ca3 < 0):
            gamma_ca3 = 0
        gamma_tot = gamma + gamma_cq + gamma_ca2 + gamma_ca3
        prop2 = 1/((sqrts - self.m_delta)**2 + gamma_tot**2)
        # Eq 2.24
        abs_pwave = self.ax2*pcm2/ppi.p()*rho*(gamma_ca2+gamma_ca3)*prop2
        # Convert to fm^-1
        abs_pwave *= HBARC**2

        # Quasielastic
        qel_pwave = self.ax2*pcm2/ppi.p()*rho*gamma*prop2
        # Convert to fm^-1
        qel_pwave *= HBARC**2
        # Eq. 3.10 mu_cm+1/2
        pauli_swave = 1 - (p*pcm-edelta*en+ef*sqrts)/(2*p*pcm)
        if(pauli_swave < 0):
            pauli_swave = 0
        elif(pauli_swave > 1):
            pauli_swave = 1
        wswave = sqrts
        if(sqrts > 6.76*HBARC):
            wswave = 6.76*HBARC
        # Eq 3.7 and 3.8    
        xi = (wswave - self.m_n - self.m_pi)/self.m_pi
        bxi = 0.21972+0.066025*xi-0.018665*xi**2
        dxi = -0.03130+0.370625*xi-0.08229*xi**2
        axi = (1+dxi)/2
        cxi = (1-dxi)/2
        sigma = 0.19753+0.06899*xi-0.01334*xi*xi
        # convert sigma to MeV^2
        sigma = sigma/self.m_pi**2
        # Eq. 3.1 
        qel_swave = rho*sigma*pauli_swave
        # Convert to fm^-1
        qel_swave *= HBARC**2

        quasielastic = (qel_pwave + qel_swave)*2*np.pi*kf**2
        absorption = (abs_pwave + abs_swave)*2*np.pi*kf**2 #+ abs_swave
        #absorption*=(self.A-1)*(self.A-2)/(self.A)**2
        #quasielastic*=(self.A-1)/(self.A)
        if(flag==1):
            res=absorption
        else:
            res=quasielastic

        return res #absorption#, quasielastic

    def density(self, radius):
        return self.rho0/(1+np.exp((radius-3.971)/0.5935))


def fit(x, a, b, c):
    return c/(1+np.exp((x-a)/b))*(HBARC/m_pi)**3


def test(k):
    return 4 * k*k/(2*np.pi)**3 * (4*np.pi)


if __name__ == '__main__':
    import pandas as pd
    from scipy.optimize import curve_fit
    from scipy.integrate import quad, dblquad,nquad

    result = quad(test, 0, 225)
    print((3*np.pi**2/2*result[0])**(1/3))

    m_pi = physics.ParticleInfo(211).mass()
    m_n = physics.ParticleInfo(2212).mass()
    tpi = 165
    epi = tpi+m_pi
    ppi = physics.Vector4(0, 0, np.sqrt(epi**2-m_pi**2), epi)

    oset = Oset(26, 30)
    #oset = Oset(20, 20)

    

    radii = np.linspace(0, 7, 20)
    rho = oset.density(radii)*(HBARC/m_pi)**3
    scatter = []
    absorbed = []
    #cos_piN=-0.95
    for radius in radii:
        rhop=oset.density(radius)
        kf=(3*np.pi**2*rhop)**(1/3)*HBARC
        absorption, err_abs = nquad(oset,[[-1.0,1.0],[0.0,kf]], args=(radius, ppi,1))
        print(absorption*3/4/np.pi/kf**3,radius,kf)  
        interact, err_int = nquad(oset,[[-1.0,1.0],[0.0,kf]], args=(radius, ppi,2))
        print(interact*3/4/np.pi/kf**3,radius,kf) 
        absorbed.append(absorption*3/4/np.pi/kf**3)
        scatter.append(interact*3/4/np.pi/kf**3)

    oset_data = pd.read_csv('oset.csv',
                            names=['Qx', 'Qy', 'Ax', 'Ay', 'rhox', 'rhoy'])

    popt, pcov = curve_fit(fit, oset_data['rhox'][:46], oset_data['rhoy'][:46])

    # Plot Oset data
    plt.plot(radii, scatter, color='tab:red')
    plt.plot(radii, absorbed, color='tab:blue')
    plt.plot(radii, rho, color='tab:green')

    # Plot our calculation
    plt.plot(oset_data['Qx'], oset_data['Qy'], color='tab:red', ls='--')
    plt.plot(oset_data['Ax'], oset_data['Ay'], color='tab:blue', ls='--')
    plt.plot(oset_data['rhox'], oset_data['rhoy'], color='tab:green', ls='--')

    plt.xlim([0,7])
    #plt.ylim([0,0.8])
    plt.show()
