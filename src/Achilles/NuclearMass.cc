// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/NuclearMass.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Exception.hh"
#include "Achilles/System.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>

#include "fmt/core.h"
#include "spdlog/spdlog.h"

namespace {

// Least-squares fit of the five Bethe-Weizsacker terms to AME2020 over A = 4-40.
// RMS 3.0 MeV on the binding energy for 8 <= A <= 40, versus ~9 MeV for the usual
// textbook coefficients; degrades above A ~ 100, which Achilles does not reach.
constexpr double aVolume = 14.039;
constexpr double aSurface = 13.921;
constexpr double aCoulomb = 0.5813;
constexpr double aAsymmetry = 17.372;
constexpr double aPairing = 8.042;

// Measured ground-state masses, read once on first use. The Bethe-Weizsacker fit
// above carries a ~3 MeV RMS error, which is the same size as the excitation
// energies it is used to test, so prefer measured values wherever the evaluation
// covers the nuclide and keep the fit only for the gaps.
constexpr const char *massTableFile = "data/ame2020_nuclear_masses.txt";

const std::map<std::pair<int, int>, double> &MassTable() {
    static const std::map<std::pair<int, int>, double> table = [] {
        std::map<std::pair<int, int>, double> masses;

        std::string path;
        try {
            path = achilles::Filesystem::FindFile(massTableFile, "NuclearMass");
        } catch(const std::exception &) {
            spdlog::warn("NuclearMass: could not find {}, falling back to the "
                         "Bethe-Weizsacker fit for every nuclide",
                         massTableFile);
            return masses;
        }

        std::ifstream data(path);
        std::string line;
        while(std::getline(data, line)) {
            if(line.empty() || line[0] == '#') continue;
            std::istringstream parser(line);
            int Z{}, A{};
            double mass{};
            if(!(parser >> Z >> A >> mass)) continue;
            masses.emplace(std::make_pair(Z, A), mass);
        }

        spdlog::debug("NuclearMass: loaded {} measured masses from {}", masses.size(), path);
        return masses;
    }();

    return table;
}

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

    // Prefer the measured mass. The fit is a stand-in for nuclides the evaluation
    // does not cover, and its ~3 MeV error is the scale of the excitation energies
    // these masses are used to compute, so the difference is not cosmetic.
    const auto &table = MassTable();
    if(const auto it = table.find({Z, A}); it != table.end()) return it->second;

    static std::once_flag warned;
    std::call_once(warned, [&] {
        spdlog::warn("NuclearMass: {} is absent from {}, falling back to the "
                     "Bethe-Weizsacker fit (reported once)",
                     achilles::NuclearName(Z, A), massTableFile);
    });

    return Z * Constant::mp + (A - Z) * Constant::mn - BindingEnergy(Z, A);
}

double achilles::SeparationEnergy(int Z, int A, bool is_proton) {
    const int dZ = is_proton ? 1 : 0;
    const double nucleon = is_proton ? Constant::mp : Constant::mn;
    return std::max(NuclearMass(Z - dZ, A - 1) + nucleon - NuclearMass(Z, A), 0.0);
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
