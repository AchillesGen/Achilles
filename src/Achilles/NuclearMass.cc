// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/NuclearMass.hh"
#include "Achilles/Constants.hh"

#include <algorithm>
#include <array>
#include <cmath>

#include "fmt/core.h"

namespace {

// Least-squares fit of the five Bethe-Weizsacker terms to AME2020 over A = 4-40.
// RMS 3.0 MeV on the binding energy for 8 <= A <= 40, versus ~9 MeV for the usual
// textbook coefficients; degrades above A ~ 100, which Achilles does not reach.
constexpr double aVolume = 14.039;
constexpr double aSurface = 13.921;
constexpr double aCoulomb = 0.5813;
constexpr double aAsymmetry = 17.372;
constexpr double aPairing = 8.042;

constexpr std::array<const char *, 119> symbols{
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf",
    "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

} // namespace

double achilles::BindingEnergy(int Z, int A) {
    // The formula only describes bound many-body systems. Everything else (free nucleons,
    // all-neutron or all-proton systems, and the A < 4 few-body states the fit excludes)
    // gets zero binding, i.e. the free-nucleon mass sum.
    if(A < 4 || Z <= 0 || Z >= A) return 0.0;

    const auto a = static_cast<double>(A);
    const auto z = static_cast<double>(Z);

    double pairing = 0.0;
    if(A % 2 == 0) pairing = (Z % 2 == 0 ? aPairing : -aPairing) / std::sqrt(a);

    const double binding = aVolume * a - aSurface * std::cbrt(a * a) -
                           aCoulomb * z * (z - 1) / std::cbrt(a) -
                           aAsymmetry * (a - 2 * z) * (a - 2 * z) / a + pairing;

    // Never let a pathological (Z, A) push the mass above the free-nucleon sum.
    return std::max(binding, 0.0);
}

double achilles::NuclearMass(int Z, int A) {
    if(A <= 0) return 0.0;
    if(A == 1) return Z == 1 ? Constant::mp : Constant::mn;

    return Z * Constant::mp + (A - Z) * Constant::mn - BindingEnergy(Z, A);
}

std::string achilles::ElementSymbol(int Z) {
    if(Z < 0 || static_cast<size_t>(Z) >= symbols.size()) return fmt::format("Z{}", Z);
    return symbols[static_cast<size_t>(Z)];
}

std::string achilles::NuclearName(int Z, int A, int L, int I) {
    std::string name = fmt::format("{}{}", ElementSymbol(Z), A);
    if(L > 0) name = std::string(static_cast<size_t>(L), 'L') + name;
    if(I > 0) name += "*";
    return name;
}
