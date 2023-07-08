#include "Achilles/Constants.hh"

#include <iostream>

// Fortran constants interface
extern "C" {

struct fconstants {
    // Constants
    double c = achilles::Constant::C;
    double hbarc = achilles::Constant::HBARC;
    double hbarc2 = achilles::Constant::HBARC2;
    double pi = M_PI;

    // Masses
    double mp = achilles::Constant::mp;
    double mn = achilles::Constant::mn;
    double mqe = achilles::Constant::mN;
    double mpi0 = achilles::Constant::mpi0;
    double mpip = achilles::Constant::mpip;
    double meta = achilles::Constant::meta;
};

void init_ (fconstants &data) {
    data.pi = M_PI;	
    data.c = achilles::Constant::C;
    data.hbarc = achilles::Constant::HBARC;
    data.hbarc2 = achilles::Constant::HBARC2;
    data.mp = achilles::Constant::mp;
    data.mn = achilles::Constant::mn;
    data.mqe = achilles::Constant::mN;
    data.mpi0 = achilles::Constant::mpi0;
    data.mpip = achilles::Constant::mpip;
    data.meta = achilles::Constant::meta;

}

}
