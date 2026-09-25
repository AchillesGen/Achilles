// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Enforcing units at the file boundary.
//
// Rule: a bare number from a file becomes a Quantity in exactly ONE place --
// here -- and only by naming the unit. A config writes the unit either as a
// plain literal, `Energy: 30 GeV`, or as `{ value: 30, unit: GeV }`; both go
// through the same check. Three enforcement levels:
//   (A) file DECLARES its unit  -> parse the string, check its dimension, apply.
//   (B) format FIXES the unit by spec (LHE=GeV, ...) -> one cited constant.
//   (C) range backstop          -> catches a 1000x MeV/GeV swap that slips (A)/(B).
// Missing/unknown/wrong-dimension units THROW; they never default silently.

#ifndef UNITS_IO_HH
#define UNITS_IO_HH

#include "Achilles/PhysicalUnits.hh"

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles::units::io {

// --- Dimension-specific unit registries (string -> Unit of that dimension) ---
// A runtime string can't carry a compile-time dimension, so each table is typed
// and a lookup that finds the wrong dimension simply isn't in the table.
inline const std::unordered_map<std::string, Unit<1>> &energy_units() {
    static const std::unordered_map<std::string, Unit<1>> t{
        {"eV", eV}, {"keV", keV}, {"MeV", MeV}, {"GeV", GeV}, {"TeV", TeV}};
    return t;
}
inline const std::unordered_map<std::string, Unit<-1>> &length_units() {
    static const std::unordered_map<std::string, Unit<-1>> t{
        {"fm", fm}, {"nm", nm}, {"m", m}, {"MeV^-1", iMeV}, {"GeV^-1", iGeV}};
    return t;
}
inline const std::unordered_map<std::string, Unit<-2>> &xsec_units() {
    static const std::unordered_map<std::string, Unit<-2>> t{
        {"MeV^-2", iMeV2}, {"GeV^-2", iGeV2}, {"fm^2", fm2}, {"b", b},  {"mb", mb},
        {"ub", ub},        {"nb", nb},        {"pb", pb},    {"fb", fb}};
    return t;
}

/// The canonical unit name for dimension D, used when writing values back out.
template <int D> constexpr const char *canonical_unit_name() {
    return D == 1 ? "MeV" : D == 2 ? "MeV^2" : D == -1 ? "MeV^-1" : D == -2 ? "MeV^-2" : "MeV^?";
}

/// Names of the units accepted for dimension D, for error messages.
template <int D> std::string known_units();

// Generic dimension-checked lookup. Wrong dimension or unknown name -> throw.
template <int D> Unit<D> unit_from_string(const std::string &name);

template <> inline std::string known_units<1>() {
    return "eV, keV, MeV, GeV, TeV";
}
template <> inline std::string known_units<-1>() {
    return "fm, nm, m, MeV^-1, GeV^-1";
}
template <> inline std::string known_units<-2>() {
    return "MeV^-2, GeV^-2, fm^2, b, mb, ub, nb, pb, fb";
}

template <> inline Unit<1> unit_from_string<1>(const std::string &name) {
    auto it = energy_units().find(name);
    if(it == energy_units().end())
        throw std::runtime_error("Units: '" + name +
                                 "' is not a known ENERGY unit (expected one of " +
                                 known_units<1>() + ")");
    return it->second;
}
template <> inline Unit<-1> unit_from_string<-1>(const std::string &name) {
    auto it = length_units().find(name);
    if(it == length_units().end())
        throw std::runtime_error("Units: '" + name +
                                 "' is not a known LENGTH unit (expected one of " +
                                 known_units<-1>() + ")");
    return it->second;
}
template <> inline Unit<-2> unit_from_string<-2>(const std::string &name) {
    auto it = xsec_units().find(name);
    if(it == xsec_units().end())
        throw std::runtime_error("Units: '" + name +
                                 "' is not a known CROSS-SECTION unit (expected one of " +
                                 known_units<-2>() + ")");
    return it->second;
}

// (A) File declares its unit as a string. Build the typed quantity, checking
//     that the declared unit actually has the dimension the caller expects.
template <int D> Quantity<D> from_declared(double value, const std::string &unit_name) {
    if(unit_name.empty())
        throw std::runtime_error("Units: missing unit, refusing to guess (expected one of " +
                                 known_units<D>() + ")");
    return value * unit_from_string<D>(unit_name); // wrong dimension -> throws
}

// (A') The same declaration written the way a physicist writes it by hand:
//      one string, "30 GeV" or "0.04fm". Whitespace between the two is
//      optional; a number with nothing after it is still a missing unit.
template <int D> Quantity<D> from_literal(const std::string &text) {
    const char *begin = text.c_str();
    char *end = nullptr;
    const double value = std::strtod(begin, &end);
    if(end == begin)
        throw std::runtime_error("Units: '" + text +
                                 "' does not start with a number (expected e.g. '30 GeV')");

    std::string unit_name(end);
    const auto first = unit_name.find_first_not_of(" \t");
    const auto last = unit_name.find_last_not_of(" \t");
    unit_name = first == std::string::npos ? "" : unit_name.substr(first, last - first + 1);

    if(unit_name.empty())
        throw std::runtime_error("Units: '" + text + "' is a bare number; name the unit, e.g. '" +
                                 text + " " + canonical_unit_name<D>() + "' (expected one of " +
                                 known_units<D>() + ")");
    return from_declared<D>(value, unit_name);
}

// (B) Format fixes the unit by specification. The unit is a compile-time arg,
//     so the dimension is checked at compile time; the citation lives in code.
template <int D> constexpr Quantity<D> from_spec(double value, Unit<D> u) noexcept {
    return value * u;
}

// (C) Physical-range backstop. Cheap tripwire for an order-of-magnitude slip.
template <int D>
Quantity<D> expect_range(Quantity<D> q, Quantity<D> lo, Quantity<D> hi, const char *what) {
    if(q < lo || q > hi)
        throw std::runtime_error(std::string("Units: value for '") + what +
                                 "' is outside its physical range (unit mistake?)");
    return q;
}

// --- Formats whose unit is fixed by their specification ---------------------

/// Les Houches Event files: all energies and momenta are in GeV.
/// (Alwall et al., "A standard format for Les Houches Event Files",
///  Comput. Phys. Commun. 176 (2007) 300, sec. 2.)
inline Energy from_lhe_energy(double v) {
    return from_spec(v, GeV);
}

/// Spectral-function and form-factor grids carry no self-description; their
/// unit comes from the data provenance and is pinned here once.
inline Energy from_spectral_grid(double v) {
    return expect_range<1>(from_spec(v, MeV), 0.0 * MeV, 5000.0 * MeV, "spectral grid energy");
}

/// Read a declared quantity and express it in `unit_name`, dispatching on
/// whichever registry knows that unit. Used where the dimension is only known
/// at run time, e.g. a validation rule that states the unit of its bounds.
inline double value_in_declared_unit(const YAML::Node &node, const std::string &unit_name);

} // namespace achilles::units::io

// ---------------------------------------------------------------------------
// YAML config. Every quantity in a config names its unit, either as a literal
//     Energy: 30 GeV
// or, equivalently, as an explicit map
//     Energy: { value: 30, unit: GeV }
// A number on its own is still refused, so a unit can never be forgotten, and
// unit_from_string<D> rejects a wrong-dimension unit at load time. Quantities
// are written back out as literals in canonical units.
// ---------------------------------------------------------------------------
namespace YAML {
template <int D> struct convert<achilles::units::Quantity<D>> {
    static bool decode(const Node &node, achilles::units::Quantity<D> &out) {
        namespace io = achilles::units::io;
        if(node.IsScalar()) {
            out = io::from_literal<D>(node.Scalar());
            return true;
        }
        if(!node.IsMap() || !node["value"] || !node["unit"])
            throw std::runtime_error("Units: expected <number> <unit> or { value: <number>, "
                                     "unit: <" +
                                     io::known_units<D>() + "> }");
        out = io::from_declared<D>(node["value"].as<double>(), node["unit"].as<std::string>());
        return true;
    }

    static Node encode(const achilles::units::Quantity<D> &q) {
        // Written as a literal in canonical units, self-describing on the way out.
        return Node(std::to_string(q.native()) + " " +
                    achilles::units::io::canonical_unit_name<D>());
    }
};
} // namespace YAML

namespace achilles::units::io {

inline double value_in_declared_unit(const YAML::Node &node, const std::string &unit_name) {
    if(energy_units().count(unit_name))
        return node.as<Quantity<1>>().in(unit_from_string<1>(unit_name));
    if(length_units().count(unit_name))
        return node.as<Quantity<-1>>().in(unit_from_string<-1>(unit_name));
    if(xsec_units().count(unit_name))
        return node.as<Quantity<-2>>().in(unit_from_string<-2>(unit_name));
    throw std::runtime_error("Units: '" + unit_name + "' is not a known unit");
}

} // namespace achilles::units::io

#endif
