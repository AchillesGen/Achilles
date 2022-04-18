#ifndef CONSTANTS_HH
#define CONSTANTS_HH

#include "Achilles/Units.hh"

namespace achilles {

namespace Constant {
    // Fundamental constants (Exact as of pdg2019)
    constexpr double C = 2.99792458e8_m / 1_s;
    constexpr double H = 6.62607015e-34_J * 1_s;
    constexpr double HBAR = H/(2*M_PI);
    constexpr double HBARC = HBAR*C;
    constexpr double HBARC2 = HBARC*HBARC*10; // mb MeV^2
    constexpr double NAVOGADRO = 6.02214076e23; // mol^-1 
    // constexpr double HBARC = 197.3269804_fm * 1_MeV;
    // constexpr double HBARC2 = 0.3893793721_mb * 1_GeV * 1_GeV;

    // Masses
    constexpr double mp = 938.27208816_MeV;
    constexpr double mn = 939.56542054_MeV;
    constexpr double mN = (mp + mn) / 2.0;
    constexpr double AMU = 931.49410248_MeV;
    constexpr double mN2 = mN*mN;
}

}

#endif
