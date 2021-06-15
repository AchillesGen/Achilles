#ifndef CONSTANTS_HH
#define CONSTANTS_HH

#include "nuchic/Units.hh"

namespace nuchic {

namespace Constant {
    // Fundamental constants
    constexpr double C = 2.9979245858E8_m / 1_s;
    constexpr double HBARC = 197.3269804_fm * 1_MeV;
    constexpr double HBARC2 = 0.3893793721_mb * 1_GeV * 1_GeV;

    // Masses
    constexpr double mp = 938.27208816_MeV;
    constexpr double mn = 939.56542052_MeV;
    constexpr double mN = (mp + mn) / 2.0;
}

}

#endif
