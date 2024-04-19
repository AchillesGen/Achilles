#ifndef CONSTANTS_HH
#define CONSTANTS_HH

// For Sherpa interface
#undef GAMMA_E

#include "Achilles/Units.hh"

namespace achilles {

namespace Constant {
// Fundamental constants (Exact as of pdg2019)
static constexpr double C = 2.99792458e8_m / 1_s;
static constexpr double H = 6.62607015e-34_J * 1_s;
static constexpr double HBAR = H / (2 * M_PI);
static constexpr double HBARC = HBAR * C;
static constexpr double HBARC2 = HBARC * HBARC * 10; // mb MeV^2
static constexpr double NAVOGADRO = 6.02214076e23;   // mol^-1
// constexpr double HBARC = 197.3269804_fm * 1_MeV;
// constexpr double HBARC2 = 0.3893793721_mb * 1_GeV * 1_GeV;
static constexpr double GAMMA_E = 0.5772156649015328606;

// Masses
static constexpr double mp = 938.27208816_MeV;
static constexpr double mn = 939.56542054_MeV;
static constexpr double mN = (mp + mn) / 2.0;
static constexpr double AMU = 931.49410248_MeV;
static constexpr double mN2 = mN * mN;
static constexpr double mpip = 139.57018_MeV;
static constexpr double mpi0 = 134.9764_MeV;
static constexpr double meta = 548.0_MeV;
static constexpr double mdelta = 1232.25_MeV;

// EW parameters
static constexpr double GF = 1.1663787e-5 / 1.0_GeV / 1.0_GeV;
static constexpr double MZ = 91.1876_GeV;
static constexpr double MW = 80.359_GeV;
static constexpr double GAMZ = 2.4952_GeV;
static constexpr double GAMW = 2.0895_GeV;
static constexpr double cos2w = MW * MW / MZ / MZ;
static constexpr double sin2w = 1 - cos2w;
static constexpr double Vud = 0.97373;
// const double alpha = sqrt(2.0)*MW*MW*GF*sin2w/M_PI;
const double alpha = 1. / 137.;
const double ee = sqrt(4 * M_PI * alpha);
const double cw = sqrt(cos2w);
const double sw = sqrt(sin2w);

} // namespace Constant

} // namespace achilles

#endif
