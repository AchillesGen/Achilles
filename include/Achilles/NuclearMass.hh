// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef NUCLEAR_MASS_HH
#define NUCLEAR_MASS_HH

#include <string>

namespace achilles {

/// Semi-empirical (Bethe-Weizsacker) binding energy in MeV.
/// Returns zero for systems the formula does not describe (A < 4, Z <= 0, Z >= A),
/// which yields an unbound free-nucleon mass from NuclearMass.
double BindingEnergy(int Z, int A);

/// Ground state nuclear (not atomic) mass in MeV for Z protons and A-Z neutrons.
double NuclearMass(int Z, int A);

/// Element symbol for a given proton number, "n" for Z == 0.
std::string ElementSymbol(int Z);

/// Name of a nuclide following the convention used in data/Particles.yml, e.g. "C11".
std::string NuclearName(int Z, int A, int L = 0, int I = 0);

} // namespace achilles

#endif
