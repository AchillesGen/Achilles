// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef CONSTANTS_HH
#define CONSTANTS_HH

// For Sherpa interface
#undef GAMMA_E

#include "Achilles/Units.hh"
#include "Achilles/UnitsFormat.hh"

#include <cmath>

namespace achilles {

namespace Constant {
// Fundamental constants (Exact as of pdg2019).
//
// These are the definitions the generator has always used; they live in
// Achilles/Units.hh so that the fm/barn unit scales derive from the very same
// hbar*c and cannot drift from these aliases.
static constexpr double C = units::kSpeedOfLight;       // fm s^-1
static constexpr double H = units::kPlanck;             // MeV s
static constexpr double HBAR = units::kHbar;            // MeV s
static constexpr double HBARC = units::kHbarC;          // MeV fm
static constexpr double HBARC2 = units::kHbarC2_MeV2mb; // mb MeV^2
static constexpr double NAVOGADRO = 6.02214076e23;      // mol^-1
static constexpr double GAMMA_E = 0.5772156649015328606;

// Masses
static constexpr units::Energy mp = 938.27208816_MeV;
static constexpr units::Energy mn = 939.56542054_MeV;
static constexpr units::Energy mN = (mp + mn) / 2.0;
static constexpr units::Energy AMU = 931.49410248_MeV;
static constexpr units::Energy2 mN2 = mN * mN;
static constexpr units::Energy mpip = 139.57018_MeV;
static constexpr units::Energy mpi0 = 134.9764_MeV;
static constexpr units::Energy meta = 548.0_MeV;
static constexpr units::Energy mdelta = 1232.25_MeV;
static constexpr units::Energy mrho = 775.8_MeV;
static constexpr units::Energy mlambda = 1115.68_MeV;
static constexpr units::Energy msigmam = 1197.45_MeV;
static constexpr units::Energy msigma0 = 1192.64_MeV;

// EW parameters
static constexpr double alpha = 1. / 137.035999084;
static constexpr units::Quantity<-2> GF = 1.1663787e-5 / 1.0_GeV / 1.0_GeV;
static constexpr double sin2w = 0.23129;
static constexpr double cos2w = 1 - sin2w;
const units::Energy MZ = units::sqrt((M_PI * alpha) / (std::sqrt(2) * GF * cos2w * sin2w));
const units::Energy MW = MZ * std::sqrt(cos2w);
static constexpr units::Energy GAMZ = 2.4952_GeV;
static constexpr units::Energy GAMW = 2.0895_GeV;
static constexpr double Vud = 0.97367;
static constexpr double Vus = 0.225;
const double ee = std::sqrt(4 * M_PI * alpha);
const double cw = std::sqrt(cos2w);
const double sw = std::sqrt(sin2w);

// Nuclear physics constants
// Nuclear densities are still carried as plain fm^-3 numbers alongside the mb
// cross sections they are used with.
static constexpr double rho0 = 0.17; // fm^-3

} // namespace Constant

} // namespace achilles

#endif
