// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "Achilles/Particle.hh"
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <utility>

namespace achilles {

class Event;

// In-medium pion absorption and quasielastic scattering follow the three papers
// of Oset et al.

// L. L. Salcedo, E. Oset, M. J. Vicente-Vacas, and
// C. Garcia-Recio, Nucl. Phys. A 484, 557 (1988)

// E. Oset and L. L. Salcedo, Nucl. Phys. A 468, 631
// (1987).

// M. J. Vicente Vacas and E. Oset, Nucl. Phys. A 568,
// 855 (1994)

class OsetPiNModel {
  public:
    OsetPiNModel() = default;

    // S + P wave
    double AbsorptionCrossSection(Event &event, std::size_t pionIndex,
                                  std::size_t nucleonIndex) const;
    double SWaveAbsorptionCrossSection(Event &event, std::size_t pionIndex,
                                       std::size_t nucleonIndex) const;
    // Keyed by {incoming pion, outgoing pion}
    std::map<std::pair<PID, PID>, double> QECrossSections(Event &event, std::size_t pionIndex,
                                                          std::size_t nucleonIndex) const;

  private:
    using Coefficients = std::array<double, 3>;
    struct CMKinematics {
        double sqrtS;
        double qcm;
    };

    static constexpr double kFermiAverage = 0.6; // <p_N^2> / kf^2

    // Delta self-energy fits: a*x^2 + b*x + c, with x = T_pi / m_pi.
    // Oset & Salcedo, Nucl. Phys. A468 (1987) 631, eq. 4.5 and Table 2; gamma = 2 beta
    static constexpr Coefficients kCQ{-5.19, 15.35, 2.06};
    static constexpr Coefficients kCA2{1.06, -6.64, 22.66};
    static constexpr Coefficients kCA3{-13.46, 46.17, -20.34};
    static constexpr Coefficients kAlpha{0.382, -1.322, 1.466};
    static constexpr Coefficients kBeta{-0.038, 0.204, 0.613};

    // piN s-wave fits in xi = (sqrt(s) - M - m_pi) / m_pi, A484 eqs. 3.7-3.8
    static constexpr Coefficients kSigma{-0.01334, 0.06889, 0.19753};
    static constexpr Coefficients kB{-0.01866, 0.06602, 0.21972};
    static constexpr Coefficients kD{-0.08229, 0.37062, -0.03130};

    static double Quadratic(double x, const Coefficients &a);
    static double CouplingFactor(double mpi);
    static CMKinematics CMFrame(const Particle &pion, double kf);
    static double RelativeSpeed(const Particle &pion, const Particle &nucleon);
    static double AbsSelfEnergy(double Tpi, double mpi, double rho);
    static double QESelfEnergy(double Tpi, double mpi, double rho);
    static double DeltaHalfWidth(double Epi, double mpi, double qcm, double kf, double sqrtS);
    static double PWaveKernel(double Epi, double mpi, double qcm, double kf, double sqrtS,
                              double rho);
    static double SWave(const Particle &pion, double rho, double vrel);
};

} // namespace achilles
