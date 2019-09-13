import numpy as np
import vegas
from nuChic.four_vector import Vec4

class PhaseSpace:
    """ Two-dimensional phase space distributions.
    
    Assumes that two beams of particles, which may possibly be monoenergetic
    interact to produce two outgoing particles.
    
    Attributes:
        beam1: first incoming beam of particles
        beam2: second incoming beam of particles
        nParticles: integer, the total number of incoming plus outgoing particles
            Default is 4.
        mass: list of floats, the masses of the particles. Defaults to massless
            particles.
    """
    
    def __init__(self,beam1,beam2,nParticles=4,mass=None):
        self._beam1 = beam1
        self._beam2 = beam2
        self._nParticles = nParticles
        if mass == None:
            self._mass = [0]*nParticles
        elif len(mass) != nParticles:
            raise Exception('Incorrect number of masses given')
        else:
            self._mass = mass

    def Generate2Body(self,x):
        """
        Generates 2-body phase space distribution.
        Args:
            x: list of two floats, "[x0, x1]", where both x0 and x1 are in the 
                unit intervals [0,1]. x0 and x1 are mapped to the angular 
                variables theta and phi, respectively
        Returns:
            (wgt, moms): the phase-space weight and the four-momenta
                wgt: float, the phase space weight
                moms: list of Vec4, the particles' 4-momenta
        """
        
        if self._nParticles != 4:
            raise Exception('More than 2 particles in the final state')

        # Generate a random value for cos(theta)
        # Use change of variables twice
        # x in [0,1] --> cos(theta) in [-1,1] --> theta in [0,pi]
        # i.e., d(cos0) = d(2x-1) = 2*dx
        dcos_theta = 2
        cos_theta = dcos_theta*x[0]-1 
        theta = np.arccos(cos_theta)

        # Generate a random value for phi
        # Use change of variables once: x in [0,1] --> phi in [0,2*pi]
        # i.e., dphi = 2*pi*dx
        dphi = 2*np.pi
        phi = dphi*x[1]

#        if self._beam1.monochromatic:
#            E1 = self._beam1.Energy()
#        else:
#            E1 = self._beam1.Energy(x[2])
#
#        if self._beam2.monochromatic:
#            E2 = self._beam2.Energy()
#        else:
#            E2 = self._beam1.Energy(x[3])

        E1 = self._beam1
        E2 = self._beam2

        # Load the mass into these variables to simplify equations
        m1 = self._mass[0]
        m2 = self._mass[1]
        m3 = self._mass[2]
        m4 = self._mass[3]

        # Find center of mass
        p1 = Vec4(E1,0,0, np.sqrt(E1**2-m1**2))
        p2 = Vec4(E2,0,0,-np.sqrt(E2**2-m2**2))
        pLab = p1 + p2
       
        # Boost to center of mass frame
        self._betaCM = pLab.BoostVector()
        p1 = p1.Boost(-self._betaCM)
        p2 = p2.Boost(-self._betaCM)

        # Calculate the momentum of the outgoing particles
        Ecm = p1.E + p2.E
        if Ecm < m3+m4:
            # Below threshold
            return 0, None
        
        # Relativistic kinematics says:
        # E_CM = E_3 + E_4
        # E_3^2 = p_CM^2 + m_3^2
        # E_4^2 = p_CM^2 + m_4^2
        # We solve these three equations for p_CM, E3, and E4
        
        # Solve for p_CM^2: p_CM^2 = (E_CM-m_3-m_4)*(E_CM+m_3-m_4)*(E_CM-m_3+m_4)*(E_CM+m_3+m_4)/(2*E_CM)^2
        # Can verify in one line with Mathematica
        pCM = (Ecm-m3-m4)*(Ecm+m3-m4)*(Ecm-m3+m4)*(Ecm+m3+m4)
        pCM = np.sqrt(pCM)/(2*Ecm)        

        # Fill the momentum
        E3 = np.sqrt(pCM**2+m3**2)
        E4 = np.sqrt(pCM**2+m4**2)
        p3 = Vec4(E3,
                  pCM*np.sin(theta)*np.cos(phi),
                  pCM*np.sin(theta)*np.sin(phi),
                  pCM*np.cos(theta))
        p4 = Vec4(E4,
                  -pCM*np.sin(theta)*np.cos(phi),
                  -pCM*np.sin(theta)*np.sin(phi),
                  -pCM*np.cos(theta))
        moms = [p1,p2,p3,p4]

        # Calculate the phase space weight
        # Assumes that incoming particles are aligned on the z-axis
        # Compare with Eq (5.32) in Schwartz's "QFT and the Standard Model"

        wgt = 1.0/(64.0*np.pi**2*Ecm**2)*p1.pz/pCM*dphi*dcos_theta

        return wgt, moms

if __name__ == '__main__':
    import vegas
    import matplotlib.pyplot as plt
    import argparse

    # Arguments to be changed on the command line
    parser = argparse.ArgumentParser()
    parser.add_argument('--E1', default=100, type=float, help='Energy of first beam in GeV')
    parser.add_argument('--E2', default=100, type=float, help='Energy of second beam in GeV')
    parser.add_argument('--m1', default=0, type=float, help='Mass of first particle in GeV')
    parser.add_argument('--m2', default=0, type=float, help='Mass of second particle in GeV')
    parser.add_argument('--m3', default=0, type=float, help='Mass of third particle in GeV')
    parser.add_argument('--m4', default=0, type=float, help='Mass of fourth particle in GeV')
    parser.add_argument('--nevents', default=10000, type=int, help='Number of events')

    # Parse input arguments
    args = vars(parser.parse_args())
    E1 = args['E1']
    E2 = args['E2']
    mass = [args['m1'],args['m2'],args['m3'],args['m4']]
    nevents = args['nevents']

    # Arrays to hold the histograms and the weight of each event
    s = []
    costheta = []
    phi = []
    wgts = []

    # Initialize the phase space class
    ps = PhaseSpace(E1,E2,4,mass)

    # hbarc2 = 3.8937966 * 10^8 pb*GeV^2 (Convert units from 1/GeV^2 to pb)
    hbarc2 = 3.8937966E8

    # Define the matrix element for e+ e- -> mu+ mu- for debugging purposes
    def Mat(moms,mass):
        """
        Matrix element for e+ e- -> mu+ mu-, intended for debugging and testing
        purposes.
        Args:
            moms: list of four-momenta as Vec4 objects. Assumes e+, e- momenta
                are the first two elements; mu+,mu- momenta are the last two
                elements 
            mass: list of the masses of the particles. Assumes e+, e- mass are
                the first two elements; mu+,mu- are mass are the last two
                elements
        Returns:
            float: the spin-averaged square of the matrix element, i.e., 
            .. math:: \frac{1}{4} \sum_{spins} |M|^2
        """        
        s = (moms[0]+moms[1]).dot(moms[0]+moms[1])
        t = (moms[0]-moms[2]).dot(moms[0]-moms[2])
        u = (moms[0]-moms[3]).dot(moms[0]-moms[3])
        
        m_e = mass[0]         # elecron mass
        m_mu = mass[2]        # muon mass
        m2 = m_e**2 + m_mu**2 # Sum of squares of masses

        alpha = 1.0/137.0     # Fine-structure constnat
        e2 = 4.0*np.pi*alpha  # Squared electric charge
        e4 = e2**2          

        # Compare to Eq (13.68) of Schwartz's "QFT and the Standard Model"
        return 2.0 * e4 * (t**2 + u**2 + 4 * s**2 * m2 - 2 * m2**2)/s**2
        

    # Generate an event given an input from VEGAS
    def GenerateEvent(x):
        """
        Generates an event using 2-body phase space
        Args:
            x: tuple of numbers (x0,x1) with x0 and x1 in the interval [0,1]
        Returns:
            wgt: float, the weighted event, i.e., the squared matrix element
                (at some point) in phase space times the phase space weight
                of the event
        Remarks: 
            depends on global variable hbarc2 to convert to Pb
            dependes on global variable fill to decide whether or not to fill
                histograms
            depends on global variables s, costheta, phi, wgts for histograms
            depends on global variable nevents to normalize the total cross section
        """
        
        wgt, event = ps.Generate2Body(x)
        if wgt == 0:
            return 0

        wgt *= Mat(event,ps._mass)*hbarc2
    
        # If on the main run, fill the histograms
        if fill:
            s.append((event[2]+event[3]).dot(event[2]+event[3]))
            costheta.append(np.cos(event[2].Theta()))
            phi.append(event[2].Phi())
    
            # Weights need to be normalized to the number of events to ensure the total cross-section is correct in the plots
            wgts.append(wgt/nevents)

        return wgt

    # Initialize VEGAS
    integ = vegas.Integrator([[0,1],[0,1]]) 

    # Preliminary run
    fill = False
    integ(GenerateEvent,nitn=10,neval=1e3)

    # Main run (fill histograms)
    fill = True
    result = integ(GenerateEvent,nitn=10,neval=nevents/10)

    # Print a summary of the VEGAS Integration results
    print("Summary of VEGAS integration results")
    print(result.summary())
    print("The best result coming from VEGAS is {0} pb".format(result))
    # Compare to the known tree-level result
    # Eq (13.78) in Schwartz's "QFT and the Standard Model" (and the discussion
    # immediately following) show that total cross section is
    # sigma = 4*pi*alpha^2 / (3* Ecm^2)
    sigma = 4.*np.pi*(1./137.)**2 / 3. / (E1 + E2)**2.
    sigma *= hbarc2 # GeV^-2 to pb
    print("The exact 1-loop answer is {0:.5f}... pb".format(sigma))
    

    # Plot the diagnostic histograms
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    ax1.hist(wgts)
    ax2.hist(costheta,weights=wgts, bins=100)
    ax3.hist(phi,weights=wgts)
    ax4.hist(s,weights=wgts)

    ax1.set_title("Weights")
    ax2.set_title("$\cos\\theta$")
    ax3.set_title("$\phi$")
    ax4.set_title("$s$")
    plt.tight_layout()
    plt.show()
