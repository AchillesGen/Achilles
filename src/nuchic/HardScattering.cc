#include "nuchic/HardScattering.hh"
#include "nuchic/Particle.hh"

nuchic::Particles nuchic::QESpectral::GeneratePhaseSpace(const std::vector<double> &rand) const {
    Particles external = {nuchic::Particle(), nuchic::Particle(),
                          nuchic::Particle(), nuchic::Particle()};

    // Incoming Electron
    FourVector mom0;
    /* Code to implement the Initial electron momentum */
    external[0].SetMomentum(mom0);

    // Incoming Nucleon
    FourVector mom1;
    /* Code to implement the Spectral Function */
    external[1].SetMomentum(mom1);

    // Outgoing Electron
    FourVector mom2;
    /* Use random numbers to setup momentum */
    external[2].SetMomentum(mom2);

    // Outgoing Nucleon
    FourVector mom3;
    /* Use random numbers to setup momentum */
    external[3].SetMomentum(mom3);

    return external;
}

double nuchic::QESpectral::CrossSection(const Particles &particles) const {
    double xsec;
    /* Calculate cross section using the external particles
     * 1. Obtain hadronic tensor
     * 2. Obtain leptonic tensor
     * 3. Contract together
     */
    return xsec;
}
