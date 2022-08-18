#include "nuchic/Constants.hh"

// Fortran constants interface
extern "C" {

struct fconstants {
    // Constants
    double c = nuchic::Constant::C;
    double hbarc = nuchic::Constant::HBARC;
    double hbarc2 = nuchic::Constant::HBARC2;

    // Masses
    double mp = nuchic::Constant::mp;
    double mn = nuchic::Constant::mn;
    double mqe = nuchic::Constant::mN;
};

void init_ (fconstants &data) {
    data.c = nuchic::Constant::C;
    data.hbarc = nuchic::Constant::HBARC;
    data.hbarc2 = nuchic::Constant::HBARC2;
    data.mp = nuchic::Constant::mp;
    data.mn = nuchic::Constant::mn;
    data.mqe = nuchic::Constant::mN;
}

}
