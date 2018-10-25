import numpy as np
import vegas
from nuChic.FourVector import Vec4

class PhaseSpace:
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
        if self._nParticles != 4:
            raise Exception('More than 2 particles in the final state')

        # Generate a random value for cos_theta
        dcos_theta = 2
        cos_theta = dcos_theta*x[0]-1
        theta = np.arccos(cos_theta)

        # Generate a random value for phi
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


        m1 = self._mass[0]
        m2 = self._mass[1]
        m3 = self._mass[2]
        m4 = self._mass[3]

        # Find center of mass
        p1 = Vec4(E1,0,0,np.sqrt(E1**2-m1**2))
        p2 = Vec4(E2,0,0,-np.sqrt(E2**2-m2**2))
        pLab = p1 + p2
       
        # Boost to center of mass frame
        betaCM = pLab.BoostVector()
        p1 = p1.Boost(pLab)
        p2 = p2.Boost(pLab)
        Ecm = p1.E + p2.E
        pCM = (Ecm-m3-m4)*(Ecm+m3-m4)*(Ecm-m3+m4)*(Ecm+m3+m4)
        pCM = np.sqrt(pCM)/(2*Ecm)        
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
        wgt = 1.0/(64.0*np.pi**2*Ecm**2)*p1.pz/pCM

        return wgt, moms

if __name__ == '__main__':
    import vegas
    import matplotlib.pyplot as plt
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--E1', default=100, type=float, help='Energy of first beam')
    parser.add_argument('--E2', default=100, type=float, help='Energy of second beam')
    parser.add_argument('--m1', default=0, type=float, help='Mass of first particle')
    parser.add_argument('--m2', default=0, type=float, help='Mass of second particle')
    parser.add_argument('--m3', default=0, type=float, help='Mass of third particle')
    parser.add_argument('--m4', default=0, type=float, help='Mass of fourth particle')
    parser.add_argument('--nevents', default=10000, type=int, help='Number of events')

    args = vars(parser.parse_args())
    E1 = args['E1']
    E2 = args['E2']
    mass = [args['m1'],args['m2'],args['m3'],args['m4']]

    s = []
    Q2Vals = []
    costheta = []
    phi = []
    wgts = []
    ps = PhaseSpace(E1,E2,4,mass)

    def Mat(moms):
        s = (moms[0]+moms[1]).dot(moms[0]+moms[1])
        t = (moms[0]-moms[2]).dot(moms[0]-moms[2])
        u = (moms[0]-moms[3]).dot(moms[0]-moms[3])

        return (t**2+u**2)/s**2

    def GenerateEvent(x):
        wgt, event = ps.Generate2Body(x)

        wgt *= Mat(event)

        if fill:
            s.append((event[2]+event[3]).dot(event[2]+event[3]))
            Q2 = (event[0]-event[2]).dot(event[0]-event[2])
            Q2Vals.append(Q2)

            costheta.append(event[2].pz/event[2].P())
            phi.append(np.arctan2(event[2].px,event[2].pz))

            wgts.append(wgt)

        return wgt

    integ = vegas.Integrator([[0,1],[0,1]]) 

    # Preliminary run
    fill = False
    integ(GenerateEvent,nitn=10,neval=1e4)

    fill = True
    result = integ(GenerateEvent,nitn=10,neval=args['nevents']/10)

    print(result.summary())

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    ax1.hist(Q2Vals,weights=wgts, bins=20)
    ax2.hist(costheta,weights=wgts, bins=20)
    ax3.hist(phi,weights=wgts)
    ax4.hist(s,weights=wgts)

    plt.tight_layout()
    plt.show()
