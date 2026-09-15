// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/CascadeInteractions/OsetPiNModel.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include <cmath>

using namespace achilles;

double OsetPiNModel::Quadratic(double x, const Coefficients &a) {
    return a[0] * x * x + a[1] * x + a[2];
}

double OsetPiNModel::CouplingFactor(double mpi) {
    return 0.36 * 4.0 * M_PI / mpi / mpi;
}

OsetPiNModel::CMKinematics OsetPiNModel::CMFrame(const Particle &pion, double kf) {
    // Average nucleon kinematics, <p_N^2> = 3/5 kf^2.
    // Simplifies calculation so one doesn't have to
    // average over the fermi sea in the calculation of the self energy
    const double pN2 = kFermiAverage * std::pow(kf, 2);

    const double pDelta2 = std::pow(pion.Momentum().P(), 2) + pN2;
    const double EDelta = pion.E() + std::sqrt(pN2 + std::pow(Constant::mN, 2));
    const double s = EDelta * EDelta - pDelta2;
    const double sqrtS = std::sqrt(s);
    const double Ecm = (s - std::pow(Constant::mN, 2) + std::pow(pion.Mass(), 2)) / (2.0 * sqrtS);
    const double qcm = std::sqrt(Ecm * Ecm - std::pow(pion.Mass(), 2));
    return {sqrtS, qcm};
}

double OsetPiNModel::RelativeSpeed(const Particle &pion, const Particle &nucleon) {
    const auto vN = nucleon.Momentum().Vec3() / nucleon.E();
    const auto vpi = pion.Momentum().Vec3() / pion.E();
    return (vpi - vN).Magnitude();
}

// A468 eq. 4.4, two- and three-body absorption terms
double OsetPiNModel::AbsSelfEnergy(double Tpi, double mpi, double rho) {
    const double x = Tpi / mpi;
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
    return abs2 + abs3;
}

// A468 eq. 4.4, quasielastic term
double OsetPiNModel::QESelfEnergy(double Tpi, double mpi, double rho) {
    const double x = Tpi / mpi;
    const double alpha = Quadratic(x, kAlpha);
    return Quadratic(x, kCQ) * std::pow(rho / Constant::rho0, alpha);
}

// A484 eqs. 2.7 and 2.13-2.14
double OsetPiNModel::DeltaHalfWidth(double Epi, double mpi, double qcm, double kf, double sqrtS) {
    const double Ef = std::sqrt(std::pow(kf, 2) + std::pow(Constant::mN, 2));
    const double EDelta =
        Epi + std::sqrt(kFermiAverage * std::pow(kf, 2) + std::pow(Constant::mN, 2));
    const double ppi = std::sqrt(Epi * Epi - mpi * mpi);
    const double pDelta = std::sqrt(std::pow(ppi, 2) + kFermiAverage * std::pow(kf, 2));
    const double ENcm = std::sqrt(qcm * qcm + std::pow(Constant::mN, 2));
    const double mu = (EDelta * ENcm - Ef * sqrtS) / pDelta / qcm;

    // Pauli-blocking reduction of the Delta decay width.
    double reduction;
    if(mu < -1.0) {
        reduction = 0.0;
    } else if(mu > 1.0) {
        reduction = 1.0;
    } else {
        reduction = (mu * mu * mu + mu + 2.0) / 4.0;
    }
    return (1.0 / 12.0 / M_PI) * CouplingFactor(mpi) * Constant::mN * qcm * qcm * qcm / sqrtS *
           reduction;
}

// A484 eqs. 2.23-2.24, common p-wave factor
double OsetPiNModel::PWaveKernel(double Epi, double mpi, double qcm, double kf, double sqrtS,
                                 double rho) {
    const double Tpi = Epi - mpi;
    const double halfWidth = DeltaHalfWidth(Epi, mpi, qcm, kf, sqrtS);
    const double re = sqrtS - Constant::mdelta;
    // Both absorption and QE self-energies enter the in-medium propagator.
    const double im = halfWidth + AbsSelfEnergy(Tpi, mpi, rho) + QESelfEnergy(Tpi, mpi, rho);
    const double prop2 = 1.0 / (re * re + im * im);

    return Constant::HBARC2 * CouplingFactor(mpi) * prop2 * std::pow(qcm, 2) / Epi;
}

double OsetPiNModel::SWave(const Particle &pion, double rho, double vrel) {
    // S-wave absorption, A484 eqs. 3.12-3.13
    const double ImB0 = 0.035;
    const double factor = 4.0 * M_PI * Constant::HBARC * 10.0 * ImB0;
    return factor / pion.E() * rho * (1.0 + pion.E() / 2.0 / Constant::mN) /
           std::pow(pion.Mass() / Constant::HBARC, 4.0) / vrel;
}

double OsetPiNModel::AbsorptionCrossSection(Event &event, std::size_t pionIndex,
                                            std::size_t nucleonIndex) const {
    const auto &pion = event.Hadrons()[pionIndex];
    const auto &nucleon = event.Hadrons()[nucleonIndex];
    const auto nucleus = event.CurrentNucleus();
    const double radius = nucleon.Position().Magnitude();
    const double rho = nucleus->ProtonRho(radius) + nucleus->NeutronRho(radius);
    const double kf = nucleus->FermiMomentum(radius, nucleon.ID());
    const auto cm = CMFrame(pion, kf);
    const double vrel = RelativeSpeed(pion, nucleon);

    // Total absorption combines the p-wave (A484 eq. 2.24) and s-wave contributions.
    const double pWave = (4.0 / 9.0) *
                         PWaveKernel(pion.E(), pion.Mass(), cm.qcm, kf, cm.sqrtS, rho) *
                         AbsSelfEnergy(pion.E() - pion.Mass(), pion.Mass(), rho) / vrel;
    return pWave + SWave(pion, rho, vrel);
}

double OsetPiNModel::SWaveAbsorptionCrossSection(Event &event, std::size_t pionIndex,
                                                 std::size_t nucleonIndex) const {
    const auto &pion = event.Hadrons()[pionIndex];
    const auto &nucleon = event.Hadrons()[nucleonIndex];
    const auto nucleus = event.CurrentNucleus();
    const double radius = nucleon.Position().Magnitude();
    const double rho = nucleus->ProtonRho(radius) + nucleus->NeutronRho(radius);
    return SWave(pion, rho, RelativeSpeed(pion, nucleon));
}

// A484 eqs. 2.16-2.19 (p-wave) and 3.3-3.8 (s-wave)
std::map<std::pair<PID, PID>, double>
OsetPiNModel::QECrossSections(Event &event, std::size_t pionIndex, std::size_t nucleonIndex) const {
    const auto &pion = event.Hadrons()[pionIndex];
    const auto &nucleon = event.Hadrons()[nucleonIndex];
    const auto nucleus = event.CurrentNucleus();
    const double radius = nucleon.Position().Magnitude();
    const double rho = nucleus->Rho(radius);
    const double kf = nucleus->FermiMomentum(radius, nucleon.ID());
    const auto cm = CMFrame(pion, kf);
    const double mpi = pion.Mass();
    const double chi = static_cast<double>(nucleus->NNeutrons() - nucleus->NProtons()) /
                       static_cast<double>(nucleus->NNucleons());

    // Q of eq. 2.19, entering the matrix of eq. 2.18 as Q/9
    const double Q = (2.0 / 3.0) * PWaveKernel(pion.E(), mpi, cm.qcm, kf, cm.sqrtS, rho) *
                     DeltaHalfWidth(pion.E(), mpi, cm.qcm, kf, cm.sqrtS);
    const double p = Q / 9.0;

    const double xi = (cm.sqrtS - Constant::mN - mpi) / mpi;
    const double s = Quadratic(xi, kSigma) / (mpi * mpi) * Constant::HBARC2;
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
