// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

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
    double mp = achilles::Constant::mp.native();
    double mn = achilles::Constant::mn.native();
    double mqe = achilles::Constant::mN.native();
    double mpi0 = achilles::Constant::mpi0.native();
    double mpip = achilles::Constant::mpip.native();
    double meta = achilles::Constant::meta.native();
    double mdelta = achilles::Constant::mdelta.native();
    double mrho = achilles::Constant::mrho.native();
};

void init_(fconstants &data) {
    data.pi = M_PI;
    data.c = achilles::Constant::C;
    data.hbarc = achilles::Constant::HBARC;
    data.hbarc2 = achilles::Constant::HBARC2;
    data.mp = achilles::Constant::mp.native();
    data.mn = achilles::Constant::mn.native();
    data.mqe = achilles::Constant::mN.native();
    data.mpi0 = achilles::Constant::mpi0.native();
    data.mpip = achilles::Constant::mpip.native();
    data.meta = achilles::Constant::meta.native();
    data.mdelta = achilles::Constant::mdelta.native();
    data.mrho = achilles::Constant::mrho.native();
}
}
