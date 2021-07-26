#ifndef CONSTANTS_HH
#define CONSTANTS_HH

#include "nuchic/Units.hh"

namespace nuchic {

namespace Constant {
    // Fundamental constants
    constexpr double C = 2.9979245858E8_m / 1_s;
    constexpr double HBARC = 200_fm * 1_MeV;
    constexpr double HBARC2 = 0.3893793721_mb * 1_GeV * 1_GeV;
    constexpr double AMU = 931.5_MeV;

    // Masses
    constexpr double mp = 938.272_MeV;
    constexpr double mn = 939.565_MeV;
    constexpr double mN = (mp + mn) / 2.0;
}

}

#endif
