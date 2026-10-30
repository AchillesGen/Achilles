// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/CascadeInteractions/OsetPiNModel.hh"
#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/PhysicalUnits.hh"
#include <cmath>

using namespace achilles;

namespace {
// Nuclear densities are carried as plain fm^-3 numbers; convert one to MeV^3
// where it meets the natural-unit cross-section formulae.
units::Quantity<3> DensityNatural(double rho) {
    const units::Length fm = 1.0_fm;
    return rho / (fm * fm * fm);
}
} // namespace

double OsetPiNModel::Quadratic(double x, const Coefficients &a) {
    return a[0] * x * x + a[1] * x + a[2];
}

units::Quantity<-2> OsetPiNModel::CouplingFactor(units::Energy mpi) {
    return 0.36 * 4.0 * M_PI / (mpi * mpi);
}

OsetPiNModel::CMKinematics OsetPiNModel::CMFrame(const Particle &pion, units::Momentum kf) {
    // Average nucleon kinematics, <p_N^2> = 3/5 kf^2.
    // Simplifies calculation so one doesn't have to
    // average over the fermi sea in the calculation of the self energy
    const units::Energy2 pN2 = kFermiAverage * kf * kf;
    const units::Energy mN = Constant::mN();
    const units::Energy mpi = pion.Mass();

    const units::Energy2 pDelta2 = pion.Momentum().P2() + pN2;
    const units::Energy EDelta = pion.E() + units::sqrt(pN2 + mN * mN);
    const units::Energy2 s = EDelta * EDelta - pDelta2;
    const units::Energy sqrtS = units::sqrt(s);
    const units::Energy Ecm = (s - mN * mN + mpi * mpi) / (2.0 * sqrtS);
    const units::Momentum qcm = units::sqrt(Ecm * Ecm - mpi * mpi);
    return {sqrtS, qcm};
}

double OsetPiNModel::RelativeSpeed(const Particle &pion, const Particle &nucleon) {
    const auto vN = nucleon.Momentum().Vec3() / nucleon.E();
    const auto vpi = pion.Momentum().Vec3() / pion.E();
    return (vpi - vN).Magnitude().native();
}

// A468 eq. 4.4, two- and three-body absorption terms
units::Energy OsetPiNModel::AbsSelfEnergy(units::Energy Tpi, units::Energy mpi, double rho) {
    const double x = (Tpi / mpi).native();
    const double r = rho / Constant::rho0;
    const double beta = Quadratic(x, kBeta);
    const double abs2 = Quadratic(x, kCA2) * std::pow(r, beta);
    double abs3 = Quadratic(x, kCA3);

    // Fit valid for 85-315 MeV; C_A3 < 0 below ~72 MeV.
    if(abs3 < 0.0) {
        abs3 = 0.0;
    } else {
        abs3 *= std::pow(r, 2.0 * beta);
    }
    return (abs2 + abs3) * units::MeV;
}

// A468 eq. 4.4, quasielastic term
units::Energy OsetPiNModel::QESelfEnergy(units::Energy Tpi, units::Energy mpi, double rho) {
    const double x = (Tpi / mpi).native();
    const double alpha = Quadratic(x, kAlpha);
    return Quadratic(x, kCQ) * std::pow(rho / Constant::rho0, alpha) * units::MeV;
}

// A484 eqs. 2.7 and 2.13-2.14
units::Energy OsetPiNModel::DeltaHalfWidth(units::Energy Epi, units::Energy mpi,
                                           units::Momentum qcm, units::Momentum kf,
                                           units::Energy sqrtS) {
    const units::Energy mN = Constant::mN();
    const units::Energy Ef = units::sqrt(kf * kf + mN * mN);
    const units::Energy EDelta = Epi + units::sqrt(kFermiAverage * kf * kf + mN * mN);
    const units::Energy2 ppi2 = Epi * Epi - mpi * mpi;
    const units::Momentum pDelta = units::sqrt(ppi2 + kFermiAverage * kf * kf);
    const units::Energy ENcm = units::sqrt(qcm * qcm + mN * mN);
    const double mu = ((EDelta * ENcm - Ef * sqrtS) / (pDelta * qcm)).native();

    // Pauli-blocking reduction of the Delta decay width.
    double reduction;
    if(mu < -1.0) {
        reduction = 0.0;
    } else if(mu > 1.0) {
        reduction = 1.0;
    } else {
        reduction = (mu * mu * mu + mu + 2.0) / 4.0;
    }
    return (1.0 / 12.0 / M_PI) * CouplingFactor(mpi) * mN * qcm * qcm * qcm / sqrtS * reduction;
}

// A484 eqs. 2.23-2.24, common p-wave factor
units::Quantity<-3> OsetPiNModel::PWaveKernel(units::Energy Epi, units::Energy mpi,
                                              units::Momentum qcm, units::Momentum kf,
                                              units::Energy sqrtS, double rho) {
    const units::Energy Tpi = Epi - mpi;
    const units::Energy halfWidth = DeltaHalfWidth(Epi, mpi, qcm, kf, sqrtS);
    const units::Energy re = sqrtS - kMdelta;
    // Both absorption and QE self-energies enter the in-medium propagator.
    const units::Energy im = halfWidth + AbsSelfEnergy(Tpi, mpi, rho) + QESelfEnergy(Tpi, mpi, rho);
    const units::Quantity<-2> prop2 = 1.0 / (re * re + im * im);

    return CouplingFactor(mpi) * prop2 * qcm * qcm / Epi;
}

units::CrossSection OsetPiNModel::SWave(const Particle &pion, double rho, double vrel) {
    // S-wave absorption, A484 eqs. 3.12-3.13. Im B0 is quoted in units of m_pi^-4.
    const double ImB0 = 0.035;
    const units::Energy mpi = pion.Mass();
    const units::Energy E = pion.E();
    const units::Quantity<-4> B0 = ImB0 / (mpi * mpi * mpi * mpi);
    return 4.0 * M_PI * B0 * DensityNatural(rho) / E *
           (1.0 + (E / (2.0 * Constant::mN())).native()) / vrel;
}

double OsetPiNModel::AbsorptionCrossSection(Event &event, std::size_t pionIndex,
                                            std::size_t nucleonIndex) const {
    const auto &pion = event.Hadrons()[pionIndex];
    const auto &nucleon = event.Hadrons()[nucleonIndex];
    const auto nucleus = event.CurrentNucleus();
    const units::Length radius = nucleon.Position().Magnitude();
    const double rho = nucleus->ProtonRho(radius) + nucleus->NeutronRho(radius);
    const units::Momentum kf = nucleus->FermiMomentum(radius, nucleon.ID());
    const auto cm = CMFrame(pion, kf);
    const double vrel = RelativeSpeed(pion, nucleon);

    // Total absorption combines the p-wave (A484 eq. 2.24) and s-wave contributions.
    const units::CrossSection pWave =
        (4.0 / 9.0) * PWaveKernel(pion.E(), pion.Mass(), cm.qcm, kf, cm.sqrtS, rho) *
        AbsSelfEnergy(pion.E() - pion.Mass(), pion.Mass(), rho) / vrel;
    return (pWave + SWave(pion, rho, vrel)).in(units::mb);
}

double OsetPiNModel::SWaveAbsorptionCrossSection(Event &event, std::size_t pionIndex,
                                                 std::size_t nucleonIndex) const {
    const auto &pion = event.Hadrons()[pionIndex];
    const auto &nucleon = event.Hadrons()[nucleonIndex];
    const auto nucleus = event.CurrentNucleus();
    const units::Length radius = nucleon.Position().Magnitude();
    const double rho = nucleus->ProtonRho(radius) + nucleus->NeutronRho(radius);
    return SWave(pion, rho, RelativeSpeed(pion, nucleon)).in(units::mb);
}

// A484 eqs. 2.16-2.19 (p-wave) and 3.3-3.8 (s-wave)
std::map<std::pair<PID, PID>, double>
OsetPiNModel::QECrossSections(Event &event, std::size_t pionIndex, std::size_t nucleonIndex) const {
    const auto &pion = event.Hadrons()[pionIndex];
    const auto &nucleon = event.Hadrons()[nucleonIndex];
    const auto nucleus = event.CurrentNucleus();
    const units::Length radius = nucleon.Position().Magnitude();
    const double rho = nucleus->Rho(radius);
    const units::Momentum kf = nucleus->FermiMomentum(radius, nucleon.ID());
    const auto cm = CMFrame(pion, kf);
    const units::Energy mpi = pion.Mass();
    // Subtract as doubles: the counts are unsigned and a proton-rich nucleus would wrap.
    const double chi =
        (static_cast<double>(nucleus->NNeutrons()) - static_cast<double>(nucleus->NProtons())) /
        static_cast<double>(nucleus->NNucleons());

    // Q of eq. 2.19, entering the matrix of eq. 2.18 as Q/9
    const units::CrossSection Q = (2.0 / 3.0) *
                                  PWaveKernel(pion.E(), mpi, cm.qcm, kf, cm.sqrtS, rho) *
                                  DeltaHalfWidth(pion.E(), mpi, cm.qcm, kf, cm.sqrtS);
    const double p = (Q / 9.0).in(units::mb);

    const double xi = ((cm.sqrtS - Constant::mN() - mpi) / mpi).native();
    const double s = (Quadratic(xi, kSigma) / (mpi * mpi)).in(units::mb);
    const double B = Quadratic(xi, kB);
    const double D = Quadratic(xi, kD);
    const double A = 0.5 * (1.0 + D);
    const double C = 0.5 * (1.0 - D);

    // Matrices of eqs. 2.18 and 3.4; pi+ <-> pi- entries vanish
    const PID pip = PID::pionp(), pi0 = PID::pion0(), pim = -PID::pionp();
    return {{{pip, pip}, s * (A - chi * B) + p * (5.0 - 4.0 * chi)},
            {{pi0, pip}, s * (1.0 - chi) * C + p * (1.0 - chi)},
            {{pip, pi0}, s * (1.0 + chi) * C + p * (1.0 + chi)},
            {{pi0, pi0}, s * D + p * 4.0},
            {{pim, pi0}, s * (1.0 - chi) * C + p * (1.0 - chi)},
            {{pi0, pim}, s * (1.0 + chi) * C + p * (1.0 + chi)},
            {{pim, pim}, s * (A + chi * B) + p * (5.0 + 4.0 * chi)}};
}
