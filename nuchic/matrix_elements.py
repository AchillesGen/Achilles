""" Basic matrix element calculations. """
from .particle import Particle
# from .pdf import PDF

MZ2 = 91.1876**2
MW2 = 80.379**2
GF = 1.166378E-5
SW2 = 0.23122
MASS_BOTTOM = 4.18
MASS_CHARM = 1.29


def main():
    """ Main function for testing. """
    from .phase_space_2d import PhaseSpace
    import vegas
    import matplotlib.pyplot as plt
    import argparse

    # Arguments to be changed on the command line
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-E',
        '--Energy',
        default=100,
        type=float,
        help='Energy of neutrino beam')
    parser.add_argument(
        '-n',
        '--nevents',
        default=10000,
        type=int,
        help='Number of events')

    # Parse input arguments
    args = vars(parser.parse_args())
    energy1 = args['Energy']
    energy2 = 1.0
    nevents = args['nevents']

    # Arrays to hold the histograms and the weight of each event
    mandel_s = []
    mandel_t = []
    mandel_u = []
    wgts = []

    # Initialize the phase space class
    phase_space = PhaseSpace(energy1, energy2, 4, [0, 1, 0, 1])

    # hbarc2 = 3.8937966 * 10^8 pb*GeV^2 * cm^2/pb (Convert units from 1/GeV^2
    # to cm^2)
    hbarc2 = 3.8937966E8 * 1E-36

    # Initialize process class for neutral processes
    process = Process(True)

    # Generate an event given an input from VEGAS
    def generate_event(x):
        wgt, event = phase_space.generate_2_body(x)
        if wgt == 0:
            return 0

        # Define the particles
        particles = [Particle(13, mom=event[0]), Particle(2212, mom=event[1]),
                     Particle(13, mom=event[2]), Particle(2212, mom=event[3])]

        # Call the Elastic scattering cross-section calculation
        wgt *= process.elastic(particles) * hbarc2

        # If on the main run, fill the histograms
        if fill:
            mandel_s.append((event[2] + event[3])**2)
            mandel_t.append((event[1] - event[3])**2)
            mandel_u.append((event[1] - event[2])**2)

            # Weights need to be normalized to the number of events to ensure
            # the total cross-section is correct in the plots
            wgts.append(wgt / nevents)

        return wgt

    # Initialize VEGAS
    integ = vegas.Integrator([[0, 1], [0, 1]])

    # Preliminary run
    fill = False
    integ(generate_event, nitn=10, neval=1e4)

    # Main run (fill histograms)
    fill = True
    result = integ(generate_event, nitn=10, neval=nevents / 10)

    # Print a summary of the VEGAS Integration results
    print(result.summary())

    # Plot the diagnostic histograms
    _, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    ax1.hist(wgts)
    ax2.hist(mandel_s, weights=wgts)
    ax3.hist(mandel_t, weights=wgts)
    ax4.hist(mandel_u, weights=wgts)
    plt.tight_layout()
    plt.show()


class Process:
    """ Base class to define an interaction process. """
    def __init__(self, neutral, pdf=None):
        self.neutral = neutral
        self.pdf = pdf

    def mott(self, particles):
        """Calculate Mott scattering matrix element."""
        raise NotImplementedError('Mott scattering is not yet implemented.')

    @staticmethod
    def elastic(particles):
        """Calculate elastic scattering matrix element from Pedro."""
        # Particle 0: incoming neutrino
        # Particle 1: incoming nucleon
        # Particle 2: outgoing neutrino
        # Particle 3: outgoing nucleon

        # charge = particles[1].charge
        # TODO fix charge to refer to a particle
        charge = 1.0
        mass_nucleon = particles[1].mass

        mandel_s = (particles[0].mom + particles[1].mom)**2
        mandel_t = (particles[0].mom - particles[2].mom)**2
        # Mandelstam u is never used?
        # u = (particles[0].mom - particles[3].mom)**2

        return GF**2 / (1 - mandel_t / MZ2)**2 * (
            32 * charge**2 * (
                mass_nucleon**4
                - 2 * mandel_s * mass_nucleon**2
                + mandel_s**2
                + 8 * mandel_t**2
                + 4 * mandel_s * mandel_t) * SW2**2
            - 8 * charge * (
                mass_nucleon**4
                - 2 * (mandel_s + 2 * mandel_t) * mass_nucleon**2
                + (mandel_s + 4 * mandel_t)**2) * SW2
            + (-mass_nucleon**2 + mandel_s + 4 * mandel_t)**2)

    def quasi_elastic(self, particles):
        """Calculate quasi-elastic matrix element."""
        raise NotImplementedError('Mott scattering is not yet implemented.')

    def dis(self, particles):
        """Calculate deep inelastic scattering (DIS) matrix element."""
        # Particle 0: incoming neutrino
        # Particle 1: incoming hadron
        # Particle 2: incoming parton
        # Particle 3: outgoing neutrino
        # Particle 4: outgoing parton

        # Calculate invariants
        # mandel_s = (particles[0].mom + particles[2].mom)**2
        mandel_t = (particles[0].mom - particles[3].mom)**2
        # mandel_u = (particles[0].mom - particles[4].mom)**2

        # Calculate the number of active fermions
        if -mandel_t > MASS_BOTTOM**2:
            nflavs = 5
        elif -mandel_t > MASS_CHARM**2:
            nflavs = 4
        else:
            nflavs = 3

        # Return 0 if parton is not allowed
        if nflavs < 5 and abs(particles[2].pid) == 5:
            return 0

        if nflavs < 4 and abs(particles[2].pid) == 4:
            return 0

        mom_transfer = particles[0].mom - particles[3].mom
        x = -mandel_t / (2 * particles[1].mom.dot(mom_transfer))
        if self.pdf is None:
            raise ValueError('PDF is not initialized.')
        f_x = self.pdf.fxQ2(particles[2].pid, x, -mandel_t)
        return f_x

    def res(self, particles):
        """Calculate resonance matrix element."""
        raise NotImplementedError('Mott scattering is not yet implemented.')


# This line will find all the processes, and store them in a dictionary,
# as long as the function name does not have a dunder in its name
PROCESSES = {k: v for (k, v) in Process.__dict__.items() if '__' not in k}

if __name__ == '__main__':
    main()
