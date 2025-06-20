#include "Achilles/OsetCrossSections.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

using achilles::OsetCrossSection;
using namespace achilles;

OsetCrossSection::OsetCrossSection() {}

double OsetCrossSection::AbsCrossSection(Event &event, size_t part1, size_t part2) const {
    // Get pion and nucleus information
    const auto &pion = event.Hadrons()[part1];
    const auto &nucleon = event.Hadrons()[part2];
    auto proton_density = event.CurrentNucleus()->ProtonRho(nucleon.Position().Magnitude());
    auto neutron_density = event.CurrentNucleus()->NeutronRho(nucleon.Position().Magnitude());
    auto total_density = proton_density + neutron_density;
    auto fermi_momentum = event.CurrentNucleus()->FermiMomentum(nucleon.Position().Magnitude());
    // spdlog::debug("Density = {}, PionE = {}, fermimom = {}", density, pion.E(), fermi_momentum);

    auto pion_kinetic_energy = pion.E() - pion.Mass();
    auto pion_momentum = pion.Momentum().P();

    // Actual nucleon kinematics included
    // auto delta_fourvector = pion.Momentum() + nucleon.Momentum();
    // auto s = delta_fourvector.M2();
    // auto sqrts = sqrt(s);

    // Let's consider an average nucleon
    double deltamomsq = pow(pion_momentum, 2) + 0.6 * pow(fermi_momentum, 2);
    double deltaE = pion.E() + sqrt(0.6 * pow(fermi_momentum, 2) + pow(Constant::mN, 2));
    double s = deltaE * deltaE - deltamomsq;
    double sqrts = sqrt(s);

    // assumes nucleon at rest?
    /*auto sqrts =
        sqrt(pion.Mass() * pion.Mass() + 2.0 * Constant::mN * pion.E() + pow(Constant::mN, 2));
    auto cms_mom = pion_momentum * Constant::mN / sqrts;
    */

    double wcm = (s - pow(Constant::mN, 2) + pow(pion.Mass(), 2)) / (2. * sqrts);
    double cms_mom2 = wcm * wcm - pow(pion.Mass(), 2);
    double cms_mom = sqrt(cms_mom2);

    // ----- ABSORPTION ----- //
    auto pXsecCommon = PXSecCommon(nucleon.E(), pion.E(), pion.Mass(), cms_mom, fermi_momentum,
                                   sqrts, total_density);

    auto vnuc = nucleon.Momentum().Vec3()/nucleon.E();
    auto vpi = pion.Momentum().Vec3()/pion.E();
    auto vrel = (vpi-vnuc).Magnitude();

    // constant factor for p-wave absorption cross section
    static const double pAbsorptionFactor = 4.0 / 9.0;

    auto self_energy_absorption = SelfEnergyAbsNN(pion_kinetic_energy, pion.Mass(), total_density) +
                                  SelfEnergyAbsNNN(pion_kinetic_energy, pion.Mass(), total_density);

    // absorption p-wave cross section (see eq. 2.24 and eq. 2.21)
    const double pXsecAbsorption = pAbsorptionFactor * pXsecCommon * self_energy_absorption / vrel;

    // constant factor for s-wave absorption cross section
    static const double sAbsorptionFactor = 4.0 * M_PI * Constant::HBARC * 10.0 * ImB0;


    // absorption s-wave cross section (see sec. 3.3)
    const double sXsecAbsorption = sAbsorptionFactor / pion.E() * total_density *
                                   (1.0 + pion.E() / 2.0 / Constant::mN) /
                                   pow(pion.Mass() / Constant::HBARC, 4.0) / vrel;

    // total absorption cross section coming from both s- and p-waves
    return pXsecAbsorption + sXsecAbsorption;
}

std::map<std::pair<PID, PID>, double> OsetCrossSection::QECrossSection(Event &event, size_t part1,
                                                                       size_t part2) const {
    // Get pion and nucleus information
    const auto &pion = event.Hadrons()[part1];
    const auto &nucleon = event.Hadrons()[part2];
    auto density = event.CurrentNucleus()->Rho(nucleon.Position().Magnitude());
    auto fermi_momentum = event.CurrentNucleus()->FermiMomentum(nucleon.Position().Magnitude());
    auto pion_momentum = pion.Momentum().P();

    double protfrac = static_cast<double>(event.CurrentNucleus()->NNeutrons() -
                                          event.CurrentNucleus()->NProtons()) /
                      static_cast<double>(event.CurrentNucleus()->NNucleons());

    // Let's consider an average nucleon with <P^2> = (3/5)kf^2
    double deltamomsq = pow(pion_momentum, 2) + 0.6 * pow(fermi_momentum, 2);
    double deltaE = pion.E() + sqrt(0.6 * pow(fermi_momentum, 2) + pow(Constant::mN, 2));
    double s = deltaE * deltaE - deltamomsq;
    double sqrts = sqrt(s);

    double wcm = (s - pow(Constant::mN, 2) + pow(pion.Mass(), 2)) / (2. * sqrts);
    double cms_mom2 = wcm * wcm - pow(pion.Mass(), 2);
    double cms_mom = sqrt(cms_mom2);

    // Compute cms quantities
    // assumes nucleon at rest?
    /*
    auto sqrts =
        sqrt(pion.Mass() * pion.Mass() + 2.0 * Constant::mN * pion.E() + pow(Constant::mN, 2));
    auto cms_mom = pion_momentum * Constant::mN / sqrts;
    */

    // constant factor for p-wave absorption cross section
    static const double pAbsorptionFactor = 4.0 / 9.0;

    double pXsecTotalQel =
        pAbsorptionFactor *
        PXSecCommon(nucleon.E(), pion.E(), pion.Mass(), cms_mom, fermi_momentum, sqrts, density) *
        ReducedHalfWidth(nucleon.E(), pion.E(), pion.Mass(), cms_mom, fermi_momentum, sqrts);

    const double ksi = (sqrts - Constant::mN - pion.Mass()) / pion.Mass();

    const double sXsecTotalQel =
        QuadraticFunction(ksi, fCoefSigma) / pow(pion.Mass(), 2) * fNormFactor;

    const double B = QuadraticFunction(ksi, fCoefB);
    const double D = QuadraticFunction(ksi, fCoefD);
    const double A = 0.5 + 0.5 * D;
    const double C = 1.0 - A; // Also (1/2(1-D))

    // Eq 3.3 + Eq. 2.18(except the paper's Qp = 1/6 * pXsecTotalQel)
    std::map<std::pair<PID, PID>, double> QelCrossSections;
    QelCrossSections[{PID::pionp(), PID::pionp()}] =
        sXsecTotalQel * (A - protfrac * B) + pXsecTotalQel * (5. - 4 * protfrac) / 6.;
    QelCrossSections[{PID::pionp(), PID::pion0()}] =
        sXsecTotalQel * (1. + protfrac) * C + pXsecTotalQel * (1 + protfrac) / 6.;
    QelCrossSections[{PID::pionp(), -PID::pionp()}] = sXsecTotalQel * 0. + 0.;
    QelCrossSections[{PID::pion0(), PID::pionp()}] =
        sXsecTotalQel * (1. - protfrac) * C + pXsecTotalQel * (1. - protfrac) / 6;
    QelCrossSections[{PID::pion0(), PID::pion0()}] = sXsecTotalQel * D + pXsecTotalQel * 4. / 6.;
    QelCrossSections[{PID::pion0(), -PID::pionp()}] =
        sXsecTotalQel * (1. + protfrac) * C + pXsecTotalQel * (1. + protfrac) / 6;
    QelCrossSections[{-PID::pionp(), PID::pionp()}] = sXsecTotalQel * 0. + pXsecTotalQel * 0.;
    QelCrossSections[{-PID::pionp(), PID::pion0()}] =
        sXsecTotalQel * (1. - protfrac) * C + pXsecTotalQel * (1. - protfrac) / 6;
    QelCrossSections[{-PID::pionp(), -PID::pionp()}] =
        sXsecTotalQel * (A + protfrac * B) + pXsecTotalQel * (5. + 4 * protfrac) / 6.;
    return QelCrossSections;
}

double OsetCrossSection::SelfEnergyAbsNN(const double pion_KE, const double pion_mass,
                                         const double density) const {
    auto x = KEFraction(pion_KE, pion_mass);
    auto beta = QuadraticFunction(x, fCoefBeta);
    auto absNN = QuadraticFunction(x, fCoefCA2);
    absNN *= pow(DensityFraction(density), beta);

    return absNN;
}

double OsetCrossSection::SelfEnergyAbsNNN(const double pion_KE, const double pion_mass,
                                          const double density) const {
    auto x = KEFraction(pion_KE, pion_mass);
    auto beta = QuadraticFunction(x, fCoefBeta);
    auto absNNN = QuadraticFunction(x, fCoefCA3);

    if(absNNN < 0.0)
        absNNN = 0.0; // it may happen for Tk < 50 MeV
    else
        absNNN *= pow(DensityFraction(density), 2.0 * beta);

    return absNNN;
}

double OsetCrossSection::SelfEnergyQE(const double pion_KE, const double pion_mass,
                                      const double density) const {
    auto x = KEFraction(pion_KE, pion_mass);
    auto alpha = QuadraticFunction(x, fCoefAlpha);
    return QuadraticFunction(x, fCoefCQ) * pow(DensityFraction(density), alpha);
}

double OsetCrossSection::PXSecCommon(const double nucE, const double pionE, const double pion_mass,
                                     const double cms_mom, const double fermimom,
                                     const double sqrts, const double density) const {
    // Compute delta prop sq
    auto pionKE = pionE - pion_mass;
    //auto pion_momentum = sqrt(pionE * pionE - pion_mass * pion_mass);
    auto reduced_halfwidth = ReducedHalfWidth(nucE, pionE, pion_mass, cms_mom, fermimom, sqrts);
    // real and imaginary part of delta denominator (see eq. 2.23)
    const double re_delta = sqrts - Constant::mdelta;
    const double im_delta = reduced_halfwidth + SelfEnergyAbsNN(pionKE, pion_mass, density) +
                            SelfEnergyAbsNNN(pionKE, pion_mass, density) +
                            SelfEnergyQE(pionKE, pion_mass, density);

    // Compute delta proopagator squared
    auto delta_propagator2 = 1.0 / (re_delta * re_delta + im_delta * im_delta);

    auto pXsecCommon = fNormFactor * CouplingFactor(pion_mass) * delta_propagator2 *
                       pow(cms_mom, 2) / pionE;

    return pXsecCommon;
}

double OsetCrossSection::ReducedHalfWidth(const double, const double pionE, const double pion_mass,
                                          const double cms_mom, const double fermimom,
                                          const double sqrts) const {
    auto fermi_energy = sqrt(pow(fermimom, 2) + pow(Constant::mN, 2));

    auto coupling_factor = CouplingFactor(pion_mass);

    // Use average nucleon momentumsq
    auto deltaE = pionE + sqrt(0.6 * pow(fermimom, 2) + pow(Constant::mN, 2));

    // real nucleon kinematics
    // auto deltaE = pionE + nucE;

    // nucleon at rest
    // auto deltaEnergy = pionE + Constant::mN;

    auto pion_mom = sqrt(pionE * pionE - pion_mass * pion_mass);
    auto deltamom = sqrt(pow(pion_mom, 2) + 0.6 * pow(fermimom, 2));
    auto n_energy_cms = sqrt(cms_mom * cms_mom + pow(Constant::mN, 2));
    auto mu0 = (deltaE * n_energy_cms - fermi_energy * sqrts) / deltamom / cms_mom;

    double delta_reduction;
    if(mu0 < -1.0)
        delta_reduction = 0.0;
    else if(mu0 > 1.0)
        delta_reduction = 1.0;
    else
        delta_reduction = (mu0 * mu0 * mu0 + mu0 + 2.) / 4.0;

    auto reduced_halfwidth = fConstFactor * coupling_factor * Constant::mN * cms_mom * cms_mom *
                             cms_mom / sqrts * delta_reduction;

    return reduced_halfwidth;
}

double avgrelvelocity(const double pion_mom, const double pionE, const double kf,
                      const double rho) {
    auto mpi = sqrt(pionE * pionE - pion_mom * pion_mom);
    auto c = pion_mom / pionE;

    auto integral = (2. / 3.) * c * (-6 + c * c) * pow(mpi, 3) +
                    (2. / 3.) * c * sqrt(kf * kf + pow(Constant::mN, 2)) *
                        ((3 + c * c) * kf * kf + (-6 + c * c) * pow(Constant::mN, 2));

    return (pow(Constant::HBARC, 3) / rho) * 4. / pow((2 * M_PI), 3) * 2 * M_PI / 3. * pow(c, -1) *
           integral;
}
