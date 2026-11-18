// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/ParticleInfo.hh"
#include "Achilles/PhysicalUnits.hh"

#include <iostream>

namespace {

// Meson and Delta masses used by the Fortran models (DCC in particular) are
// tuned against these values, so they are pinned here rather than read from
// Particles.yml, which users are free to retune.
constexpr double kMpi0 = 134.9764;  // MeV
constexpr double kMpip = 139.57018; // MeV
constexpr double kMeta = 548.0;     // MeV
constexpr double kMdelta = 1232.25; // MeV
constexpr double kMrho = 775.8;     // MeV

} // namespace

// Fortran constants interface
extern "C" {

struct fconstants {
    // Constants
    double c = achilles::Constant::C;
    double hbarc = achilles::Constant::HBARC;
    double hbarc2 = achilles::Constant::HBARC2;
    double pi = M_PI;

    // Masses
    double mp = achilles::Constant::mp().native();
    double mn = achilles::Constant::mn().native();
    double mqe = achilles::Constant::mN().native();
    double mpi0 = kMpi0;
    double mpip = kMpip;
    double meta = kMeta;
    double mdelta = kMdelta;
    double mrho = kMrho;
};

void init_(fconstants &data) {
    data.pi = M_PI;
    data.c = achilles::Constant::C;
    data.hbarc = achilles::Constant::HBARC;
    data.hbarc2 = achilles::Constant::HBARC2;
    data.mp = achilles::Constant::mp().native();
    data.mn = achilles::Constant::mn().native();
    data.mqe = achilles::Constant::mN().native();
    data.mpi0 = kMpi0;
    data.mpip = kMpip;
    data.meta = kMeta;
    data.mdelta = kMdelta;
    data.mrho = kMrho;
}
}
