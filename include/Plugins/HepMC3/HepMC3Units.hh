// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// HepMC3 files declare their own units per event. Read them; never assume.

#ifndef HEPMC3_UNITS_HH
#define HEPMC3_UNITS_HH

#include "Achilles/FourVector.hh"
#include "Achilles/UnitsIO.hh"

#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"

namespace achilles::units::io {

/// The energy unit this event declares.
inline Unit<1> hepmc_energy_unit(const HepMC3::GenEvent &e) {
    switch(e.momentum_unit()) {
    case HepMC3::Units::MEV:
        return MeV;
    case HepMC3::Units::GEV:
        return GeV;
    }
    throw std::runtime_error("HepMC3: event declares an unknown momentum unit");
}

/// The length unit this event declares.
inline Unit<-1> hepmc_length_unit(const HepMC3::GenEvent &e) {
    // 1 mm = 1e12 fm, 1 cm = 1e13 fm.
    switch(e.length_unit()) {
    case HepMC3::Units::MM:
        return Unit<-1>{fm.scale() * 1.0e12, "mm"};
    case HepMC3::Units::CM:
        return Unit<-1>{fm.scale() * 1.0e13, "cm"};
    }
    throw std::runtime_error("HepMC3: event declares an unknown length unit");
}

/// Read an incoming HepMC3 momentum into a typed Achilles FourVector.
inline achilles::FourVector momentum_from_hepmc(const HepMC3::FourVector &p,
                                                const HepMC3::GenEvent &e) {
    const auto u = hepmc_energy_unit(e);
    return {p.e() * u, p.px() * u, p.py() * u, p.pz() * u};
}

/// Read an incoming HepMC3 vertex position into a typed Achilles FourPosition.
inline achilles::FourPosition position_from_hepmc(const HepMC3::FourVector &x,
                                                  const HepMC3::GenEvent &e) {
    const auto u = hepmc_length_unit(e);
    return {x.t() * u, x.x() * u, x.y() * u, x.z() * u};
}

} // namespace achilles::units::io

#endif
