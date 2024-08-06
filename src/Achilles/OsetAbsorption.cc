#include "Achilles/OsetAbsorption.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

using achilles::OsetAbsCrossSection;

OsetAbsCrossSection::OsetAbsCrossSection() {}
double OsetAbsCrossSection::AbsorptionCrossSection(Event &event, size_t part1, size_t) const {
    // Get pion and nucleus information
    const auto &pion = event.Hadrons()[part1];
    auto density = event.CurrentNucleus()->Rho(pion.Position().Magnitude());
    auto fermi_momentum = event.CurrentNucleus()->FermiMomentum(pion.Position().Magnitude());
    auto fermi_energy = sqrt(pow(fermi_momentum, 2) + pow(Constant::mN, 2));
    auto pion_kinetic_energy = pion.E() - pion.Mass();
    auto pion_momentum = pion.Momentum().P();

    // Compute cms quantities
    // assumes nucleon at rest?
    auto sqrts =
        sqrt(pion.Mass() * pion.Mass() + 2.0 * Constant::mN * pion.E() + pow(Constant::mN, 2));
    auto momentum_cms = pion_momentum * Constant::mN / sqrts;
    auto momentum_cms2 = momentum_cms * momentum_cms;

    // Set up Delta couplings
    static const double constFactor = 1.0 / 12.0 / M_PI;
    auto coupling_factor =
        fCouplingConstant / pion.Mass() / pion.Mass(); // (f*/m_pi)^2, e.g. eq. 2.6

    // Reduce delta width due to pauli blocking see eq. 2.7 and sec. 2.3
    auto reduced_halfwidth =
        constFactor * coupling_factor * Constant::mN * momentum_cms * momentum_cms2 / sqrts *
        DeltaReduction(pion.E(), pion.Mass(), momentum_cms, fermi_energy, sqrts);

    // Calculate QEL and Abs part of deleta self energy
    const double x = pion_kinetic_energy / pion.Mass();
    const double densityFraction = density / fNormalDensity;

    const double alpha = QuadraticFunction(x, fCoefAlpha);
    const double beta = QuadraticFunction(x, fCoefBeta);

    double absNN = QuadraticFunction(x, fCoefCA2);
    double absNNN = QuadraticFunction(x, fCoefCA3);

    absNN *= pow(densityFraction, beta);

    if(absNNN < 0.0)
        absNNN = 0.0; // it may happen for Tk < 50 MeV
    else
        absNNN *= pow(densityFraction, 2.0 * beta);

    auto self_energy_absorption = absNN + absNNN;

    // this one goes to the delta propagator
    auto self_energy_total =
        absNN + absNNN + QuadraticFunction(x, fCoefCQ) * pow(densityFraction, alpha);

    // real and imaginary part of delta denominator (see eq. 2.23)
    const double re_delta = sqrts - Constant::mdelta;
    const double im_delta = reduced_halfwidth + self_energy_total;

    // Compute delta proopagator squared
    auto delta_propagator2 = 1.0 / (re_delta * re_delta + im_delta * im_delta);

    // common part for all p-wave cross sections
    const double pXsecCommon =
        fNormFactor * coupling_factor * delta_propagator2 * momentum_cms2 / pion_momentum;

    // ----- ABSORPTION ----- //

    // constant factor for p-wave absorption cross section
    static const double pAborptionFactor = 4.0 / 9.0;

    // absorption p-wave cross section (see eq. 2.24 and eq. 2.21)
    const double pXsecAbsorption = pAborptionFactor * pXsecCommon * self_energy_absorption;

    // constant factor for s-wave absorption cross section
    static const double sAbsorptionFactor = 4.0 * M_PI * Constant::HBARC * 10.0 * ImB0;

    // absorption s-wave cross section (see sec. 3.3)
    const double sXsecAbsorption = sAbsorptionFactor / pion_momentum * density *
                                   (1.0 + pion.E() / 2.0 / Constant::mN) /
                                   pow(pion.Mass() / Constant::HBARC, 4.0);

    // total absorption cross section coming from both s- and p-waves
    return pXsecAbsorption + sXsecAbsorption;
}

double OsetAbsCrossSection::DeltaReduction(const double pion_E, const double pion_mass,
                                           const double momentum_cms, const double fermi_energy,
                                           const double sqrts) const {
    // assuming nucleon at rest
    const double deltaEnergy = pion_E + Constant::mN;

    const double pion_p = sqrt(pion_E * pion_E - pion_mass * pion_mass);
    // nucleon energy in CMS
    const double energy_cms = sqrt(momentum_cms * momentum_cms + pow(Constant::mN, 2));
    // fermiEnergy calculated from density
    const double mu0 = (deltaEnergy * energy_cms - fermi_energy * sqrts) / pion_p / momentum_cms;

    // see eq. 2.13
    if(mu0 < -1.0) return 0.0;
    if(mu0 > 1.0) return 1.0;

    return (mu0 * mu0 * mu0 + mu0 + 2) / 4.0;
}
