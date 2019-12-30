""" Implement a class for generating 2-d phase space points. """

import numpy as np
import vegas
from .four_vector import Vec4
from .constants import HBARC2


def main():
    """ Function for testing basic functionality. """
    import matplotlib.pyplot as plt
    import argparse

    # Arguments to be changed on the command line
    parser = argparse.ArgumentParser()
    parser.add_argument('--E1', default=100, type=float,
                        help='Energy of first beam in GeV')
    parser.add_argument('--E2', default=100, type=float,
                        help='Energy of second beam in GeV')
    parser.add_argument(
        '--m1',
        default=0,
        type=float,
        help='Mass of first particle in GeV')
    parser.add_argument(
        '--m2',
        default=0,
        type=float,
        help='Mass of second particle in GeV')
    parser.add_argument(
        '--m3',
        default=0,
        type=float,
        help='Mass of third particle in GeV')
    parser.add_argument(
        '--m4',
        default=0,
        type=float,
        help='Mass of fourth particle in GeV')
    parser.add_argument(
        '--nevents',
        default=10000,
        type=int,
        help='Number of events')

    # Parse input arguments
    args = vars(parser.parse_args())

    # Arrays to hold the histograms and the weight of each event
    plots = {'s': [],
             'costheta': [],
             'phi': [],
             'wgts': []}

    # Initialize the phase space class
    phase_space = PhaseSpace(args['E1'], args['E2'], 4,
                             [args['m1'], args['m2'], args['m3'], args['m4']])

    # Define the matrix element for e+ e- -> mu+ mu- for debugging purposes
    def mat(moms, mass):
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
            .. math:: \\frac{1}{4} \\sum_{spins} |M|^2
        """
        mandelstam = [(moms[0] + moms[1]).dot(moms[0] + moms[1]),
                      (moms[0] - moms[2]).dot(moms[0] - moms[2]),
                      (moms[0] - moms[3]).dot(moms[0] - moms[3])]

        m_2 = mass[0]**2 + mass[2]**2  # Sum of squares of masses

        alpha = 1.0 / 137.0     # Fine-structure constnat
        e_2 = 4.0 * np.pi * alpha  # Squared electric charge
        e_4 = e_2**2

        # Compare to Eq (13.68) of Schwartz's "QFT and the Standard Model"
        return 2.0 * e_4 * (mandelstam[1]**2 + mandelstam[2]**2
                            + 4 * mandelstam[0]**2 * m_2 - 2 * m_2**2) \
            / mandelstam[0]**2

    # Generate an event given an input from VEGAS
    def generate_event(x):
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
            depends on global variable nevents to normalize the total cross
                section
        """

        wgt, event = phase_space.generate_2_body(x)
        if wgt == 0:
            return 0

        wgt *= mat(event, phase_space.mass) * HBARC2

        # If on the main run, fill the histograms
        if fill:
            plots['s'].append((event[2] + event[3]).dot(event[2] + event[3]))
            plots['costheta'].append(np.cos(event[2].Theta()))
            plots['phi'].append(event[2].Phi())

            # Weights need to be normalized to the number of events to ensure
            # the total cross-section is correct in the plots
            plots['wgts'].append(wgt / args['nevents'])

        return wgt

    # Initialize VEGAS
    integ = vegas.Integrator([[0, 1], [0, 1]])

    # Preliminary run
    fill = False
    integ(generate_event, nitn=10, neval=1e3)

    # Main run (fill histograms)
    fill = True
    result = integ(generate_event, nitn=10, neval=args['nevents'] / 10)

    # Print a summary of the VEGAS Integration results
    print("Summary of VEGAS integration results")
    print(result.summary())
    print("The best result coming from VEGAS is {0} pb".format(result))
    # Compare to the known tree-level result
    # Eq (13.78) in Schwartz's "QFT and the Standard Model" (and the discussion
    # immediately following) show that total cross section is
    # sigma = 4*pi*alpha^2 / (3* energy_cm^2)
    sigma = 4. * np.pi * (1. / 137.)**2 / 3. / (args['E1'] + args['E2'])**2.
    sigma *= HBARC2  # GeV^-2 to pb
    print("The exact 1-loop answer is {0:.5f}... pb".format(sigma))

    # Plot the diagnostic histograms
    def plot():
        _, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
        ax1.hist(plots['wgts'])
        ax2.hist(plots['costheta'], weights=plots['wgts'], bins=100)
        ax3.hist(plots['phi'], weights=plots['wgts'])
        ax4.hist(plots['s'], weights=plots['wgts'])

        ax1.set_title("Weights")
        ax2.set_title(r"$\cos\\theta$")
        ax3.set_title(r"$\phi$")
        ax4.set_title("$s$")
        plt.tight_layout()
        plt.show()

    plot()


class PhaseSpace:
    """ Two-dimensional phase space distributions.

    Assumes that two beams of particles, which may possibly be monoenergetic
    interact to produce two outgoing particles.

    Attributes:
        beam1: first incoming beam of particles
        beam2: second incoming beam of particles
        nParticles: integer, the total number of incoming plus outgoing
            particles. Default is 4.
        mass: list of floats, the masses of the particles. Defaults to massless
            particles.
    """

    def __init__(self, beam1, beam2, n_particles=4, mass=None):
        self.beam1 = beam1
        self.beam2 = beam2
        self.n_particles = n_particles
        if mass is None:
            self.mass = [0] * n_particles
        elif len(mass) != n_particles:
            raise Exception('Incorrect number of masses given')
        else:
            self.mass = mass

    @staticmethod
    def _angle_map(x):
        # Generate a random value for cos(theta)
        # Use change of variables twice
        # x in [0,1] --> cos(theta) in [-1,1] --> theta in [0,pi]
        # i.e., d(cos0) = d(2x-1) = 2*dx
        dcos_theta = 2
        cos_theta = dcos_theta * x[0] - 1
        theta = np.arccos(cos_theta)

        # Generate a random value for phi
        # Use change of variables once: x in [0,1] --> phi in [0,2*pi]
        # i.e., dphi = 2*pi*dx
        dphi = 2 * np.pi
        phi = dphi * x[1]

        return [theta, phi], dcos_theta*dphi

    def __call__(self, x):
        """ Alias for generate_2_body. """
        return self.generate_2_body(x)

    def generate_2_body(self, x):
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

        if self.n_particles != 4:
            raise Exception('More than 2 particles in the final state')

        angle, angle_wgt = self._angle_map(x)

#        if self._beam1.monochromatic:
#            energy1 = self._beam1.Energy()
#        else:
#            energy1 = self._beam1.Energy(x[2])
#
#        if self._beam2.monochromatic:
#            energy2 = self._beam2.Energy()
#        else:
#            energy2 = self._beam1.Energy(x[3])

        energy1 = self.beam1
        energy2 = self.beam2

        # Load the mass into these variables to simplify equations
        mass1 = self.mass[0]
        mass2 = self.mass[1]
        mass3 = self.mass[2]
        mass4 = self.mass[3]

        # Find center of mass
        p_1 = Vec4(energy1, 0, 0, np.sqrt(energy1**2 - mass1**2))
        p_2 = Vec4(energy2, 0, 0, -np.sqrt(energy2**2 - mass2**2))

        # Boost to center of mass frame
        _beta_cm = (p_1 + p_2).boost_vector()
        p_1 = p_1.boost(-_beta_cm)
        p_2 = p_2.boost(-_beta_cm)

        # Calculate the momentum of the outgoing particles
        energy_cm = p_1.energy + p_2.energy
        if energy_cm < mass3 + mass4:
            # Below threshold
            return 0, None

        # Relativistic kinematics says:
        # E_CM = E_3 + E_4
        # E_3^2 = p_CM^2 + m_3^2
        # E_4^2 = p_CM^2 + m_4^2
        # We solve these three equations for p_CM, E3, and E4

        # Solve for p_CM^2:
        #           (E_CM-m_3-m_4)*(E_CM+m_3-m_4)*(E_CM-m_3+m_4)*(E_CM+m_3+m_4)
        #  p_CM^2 = -----------------------------------------------------------
        #                                   (2*E_CM)^2
        # Can verify in one line with Mathematica
        p_cm = (energy_cm - mass3 - mass4) * (energy_cm + mass3 - mass4) * \
            (energy_cm - mass3 + mass4) * (energy_cm + mass3 + mass4)
        p_cm = np.sqrt(p_cm) / (2 * energy_cm)

        # Fill the momentum
        p_3 = Vec4(p_cm**2 + mass3**2,
                   p_cm * np.sin(angle[0]) * np.cos(angle[1]),
                   p_cm * np.sin(angle[0]) * np.sin(angle[1]),
                   p_cm * np.cos(angle[0]))
        p_4 = Vec4(p_cm**2 + mass4**2,
                   -p_cm * np.sin(angle[0]) * np.cos(angle[1]),
                   -p_cm * np.sin(angle[0]) * np.sin(angle[1]),
                   -p_cm * np.cos(angle[0]))
        moms = [p_1, p_2, p_3, p_4]

        # Calculate the phase space weight
        # Assumes that incoming particles are aligned on the z-axis
        # Compare with Eq (5.32) in Schwartz's "QFT and the Standard Model"

        wgt = 1.0 / (64.0 * np.pi**2 * energy_cm**2) * \
            p_1.p_z / p_cm * angle_wgt

        return wgt, moms


if __name__ == '__main__':
    main()
