// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Column-level units for homogeneous data tables (e.g. data/Particles.yml).
//
// The unit is a property of the COLUMN, not each cell. It is resolved ONCE per
// file and applied to every row -- so the compact positional layout is
// untouched and there is exactly one place a mistake can enter.
//
//   Particles:
//       units: { mass: MeV, width: MeV }   # <-- declared once (optional)
//       # [pid, mass, width, icharge, ...]
//       - Particle: [2212, 938.27, 0.0, ...]
//
// If the block is omitted, the reader falls back to Achilles' DOCUMENTED
// convention (MeV) -- that convention is itself the declaration, not a guess.
// A present block overrides it and is dimension-checked at load.

#ifndef UNITS_SCHEMA_HH
#define UNITS_SCHEMA_HH

#include "Achilles/UnitsIO.hh"

namespace achilles::units::io {

// Units for the columns of the particle table. Defaults encode the convention.
struct ParticleColumnUnits {
    Unit<1> mass{MeV};
    Unit<1> width{MeV};
};

/// Resolve the column units from an optional `units:` block. A key that is
/// present is dimension-checked; a key that is absent keeps the convention.
inline ParticleColumnUnits resolve_particle_units(const YAML::Node &units_node) {
    ParticleColumnUnits u;                           // start from the documented convention (MeV)
    if(!units_node || units_node.IsNull()) return u; // absent -> the convention stands
    if(!units_node.IsMap())
        throw std::runtime_error(
            "Particles: 'units' must be a map, e.g. { mass: MeV, width: MeV }");
    if(units_node["mass"]) u.mass = unit_from_string<1>(units_node["mass"].as<std::string>());
    if(units_node["width"]) u.width = unit_from_string<1>(units_node["width"].as<std::string>());
    return u;
}

// Apply the resolved column units to a raw row cell, with a range backstop.
inline Mass particle_mass(double raw, const ParticleColumnUnits &u) {
    return expect_range<1>(raw * u.mass, 0.0 * MeV, 1.0e6 * MeV, "particle mass");
}
inline Energy particle_width(double raw, const ParticleColumnUnits &u) {
    return expect_range<1>(raw * u.width, 0.0 * MeV, 1.0e6 * MeV, "particle width");
}

} // namespace achilles::units::io

#endif
