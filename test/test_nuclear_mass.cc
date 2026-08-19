// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/NuclearMass.hh"

using achilles::BindingEnergy;
using achilles::NuclearMass;
using Catch::Matchers::WithinAbs;

TEST_CASE("Nuclear masses reproduce AME2020 within the fit accuracy", "[NuclearMass]") {
    // AME2020 ground state nuclear masses in MeV.
    CHECK_THAT(NuclearMass(6, 12), WithinAbs(11174.86, 5.0));
    CHECK_THAT(NuclearMass(8, 16), WithinAbs(14895.08, 5.0));
    CHECK_THAT(NuclearMass(18, 40), WithinAbs(37215.53, 5.0));
    CHECK_THAT(NuclearMass(6, 11), WithinAbs(10254.02, 5.0));
    CHECK_THAT(NuclearMass(5, 11), WithinAbs(10252.55, 5.0));
}

TEST_CASE("Free nucleons bypass the mass formula", "[NuclearMass]") {
    CHECK(NuclearMass(1, 1) == achilles::Constant::mp);
    CHECK(NuclearMass(0, 1) == achilles::Constant::mn);
    CHECK(BindingEnergy(1, 1) == 0.0);
    CHECK(BindingEnergy(0, 1) == 0.0);
}

TEST_CASE("Unbound systems get the free nucleon sum", "[NuclearMass]") {
    // No protons, no Coulomb partner, or too few nucleons for the fit domain.
    CHECK(BindingEnergy(0, 2) == 0.0);
    CHECK(NuclearMass(0, 2) == 2 * achilles::Constant::mn);
    CHECK(BindingEnergy(2, 2) == 0.0);
    CHECK(NuclearMass(2, 2) == 2 * achilles::Constant::mp);
    CHECK(BindingEnergy(1, 2) == 0.0);
    CHECK(NuclearMass(0, 0) == 0.0);
}

TEST_CASE("Binding energy is never negative", "[NuclearMass]") {
    for(int A = 0; A < 300; ++A) {
        for(int Z = 0; Z <= A; ++Z) { REQUIRE(BindingEnergy(Z, A) >= 0.0); }
    }
}

TEST_CASE("Mass grows with nucleon number", "[NuclearMass]") {
    for(int Z = 1; Z < 20; ++Z) {
        for(int A = Z + 1; A < 60; ++A) { REQUIRE(NuclearMass(Z, A) > NuclearMass(Z, A - 1)); }
    }
}

TEST_CASE("Element symbols and nuclide names", "[NuclearMass]") {
    CHECK(achilles::ElementSymbol(0) == "n");
    CHECK(achilles::ElementSymbol(6) == "C");
    CHECK(achilles::ElementSymbol(18) == "Ar");
    CHECK(achilles::NuclearName(6, 11) == "C11");
    CHECK(achilles::NuclearName(6, 11, 0, 1) == "C11*");
    CHECK(achilles::NuclearName(6, 12, 1) == "LC12");
}
