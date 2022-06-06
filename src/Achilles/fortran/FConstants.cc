#include "Achilles/Constants.hh"

#include <iostream>

// Fortran constants interface
extern "C" {

struct fconstants {
    // Constants
    double c = achilles::Constant::C;
    double hbarc = achilles::Constant::HBARC;
    double hbarc2 = achilles::Constant::HBARC2;

    // Masses
    double mp = achilles::Constant::mp;
    double mn = achilles::Constant::mn;
    double mqe = achilles::Constant::mN;
};

void init_ (fconstants &data) {
    data.c = achilles::Constant::C;
    data.hbarc = achilles::Constant::HBARC;
    data.hbarc2 = achilles::Constant::HBARC2;
    data.mp = achilles::Constant::mp;
    data.mn = achilles::Constant::mn;
    data.mqe = achilles::Constant::mN;
}

}
