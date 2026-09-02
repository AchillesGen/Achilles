// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Tests for the strong-typed natural-units layer (Achilles/Units.hh).

#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/Units.hh"
#include "Achilles/UnitsFormat.hh"

#include <type_traits>

using namespace achilles::units;
using namespace achilles::units::literals;
using Catch::Approx;

// Physics anchors: constants that, if silently changed, shift every result.
TEST_CASE("Unit anchors match PDG / legacy constants", "[Units]") {
    // hbar*c is built from the SI-exact c and h, exactly as Constants.hh always has;
    // the PDG's rounded 197.3269804 agrees to 3e-10.
    CHECK(kHbarC == Approx(197.32698045930246).epsilon(1e-15));
    CHECK(kHbarC == Approx(197.3269804).epsilon(1e-9));
    CHECK(kHbarC2_MeV2mb == Approx(389379.37).epsilon(1e-8));          // old Constants::HBARC2
    CHECK((1.0 * iGeV2).in(mb) == Approx(0.3893793721).epsilon(1e-7)); // PDG
    CHECK((1.0 * fm2).in(mb) == Approx(10.0).epsilon(1e-12));          // barn definition
    CHECK((1.0 * b).in(fm2) == Approx(100.0).epsilon(1e-12));

    // Constants.hh must be a view onto the same bridge constant, not a copy.
    CHECK(achilles::Constant::HBARC == kHbarC);
    CHECK(achilles::Constant::HBARC2 == kHbarC2_MeV2mb);
    CHECK(achilles::Constant::HBARC2 == Approx(389379.37).epsilon(1e-8));
    CHECK(achilles::Constant::HBARC * achilles::Constant::HBARC * 10 ==
          Approx(achilles::Constant::HBARC2).epsilon(1e-12));
}

TEST_CASE("Unit round-trips are exact", "[Units]") {
    CHECK((2.2_GeV).in(GeV) == Approx(2.2));
    CHECK((2.2_GeV).in(MeV) == Approx(2200.0));
    CHECK((1.2_fm).in(fm) == Approx(1.2));
    CHECK((3.0e-11 * iMeV2).in(iMeV2) == Approx(3.0e-11));
    CHECK((1.0_fm).native() == Approx(1.0 / kHbarC).epsilon(1e-15));

    SECTION("Cross-section units interconvert through the same barn scale") {
        CHECK((1.0_mb).in(nb) == Approx(1.0e6));
        CHECK((2.5_nb).in(nb) == Approx(2.5));
        CHECK((1.0_mb).in(iMeV2) == Approx(1.0 / kHbarC2_MeV2mb).epsilon(1e-12));
        CHECK((1.0 * iGeV2).in(iMeV2) == Approx(1.0e-6));
    }

    SECTION("Inverse-energy units") {
        CHECK((1.0 * iGeV).in(iMeV) == Approx(1.0e-3));
        CHECK((1.0 * m).in(fm) == Approx(1.0e15));
        CHECK((1.0_m).in(fm) == Approx(1.0e15));
    }
}

TEST_CASE("Quantity algebra tracks mass dimension", "[Units]") {
    const Energy e = 2.0_GeV;
    const Length r = 1.0_fm;

    STATIC_REQUIRE(std::is_same<decltype(e * e), Energy2>::value);
    STATIC_REQUIRE(std::is_same<decltype(e * r), Dimensionless>::value);
    STATIC_REQUIRE(std::is_same<decltype(1.0 / e), Length>::value);
    STATIC_REQUIRE(std::is_same<decltype(sqrt(e * e)), Energy>::value);
    STATIC_REQUIRE(std::is_same<decltype(e / r), Energy2>::value);

    CHECK((e * e).in(GeV2) == Approx(4.0));
    CHECK(sqrt(e * e).in(GeV) == Approx(2.0));
    CHECK((e * r).native() == Approx(2000.0 / kHbarC));
    CHECK((1.0 / e).in(fm) == Approx(kHbarC / 2000.0));
    CHECK((-e).in(GeV) == Approx(-2.0));
    CHECK(abs(-e) == e);
    CHECK(e + e == 4.0_GeV);
    CHECK(e - e == Energy{});
    CHECK(e / 2.0 == 1.0_GeV);
    CHECK(2.0 * e == 4.0_GeV);
    CHECK(0.5_GeV < e);
}

TEST_CASE("sqrt is constexpr and agrees with std::sqrt at runtime", "[Units]") {
    constexpr Energy ce = sqrt(Energy2{4.0e6});
    STATIC_REQUIRE(ce.native() == 2000.0);

    Energy2 s{1234.567}; // runtime value: takes the std::sqrt path
    CHECK(sqrt(s).native() == std::sqrt(1234.567));
    CHECK(sqrt(Energy2{}).native() == 0.0);
}

TEST_CASE("Quantity is a zero-overhead double", "[Units]") {
    STATIC_REQUIRE(sizeof(Energy) == sizeof(double));
    STATIC_REQUIRE(std::is_trivially_copyable<Energy>::value);
    STATIC_REQUIRE(std::is_standard_layout<Energy>::value);
    STATIC_REQUIRE(!std::is_convertible<double, Energy>::value);
    STATIC_REQUIRE(!std::is_convertible<Energy, double>::value);
    STATIC_REQUIRE(std::is_constructible<Energy, double>::value);
}

TEST_CASE("Quantity formats in canonical units by default", "[Units]") {
    CHECK(fmt::format("{}", 2.2_GeV) == "2200 MeV");
    CHECK(fmt::format("{:.4f}", 1.0 * iMeV2) == "1.0000 MeV^-2");
    CHECK(fmt::format("{}", Dimensionless{0.5}) == "0.5");

    // A length prints as MeV^-1 because that is how it is stored.
    CHECK(fmt::format("{:.5f}", 1.2_fm) == "0.00608 MeV^-1");
    CHECK(fmt::format("{:.3f}", in(1.2_fm, fm)) == "1.200 fm");
    CHECK(fmt::format("{:.3f}", in(2.2_GeV, GeV)) == "2.200 GeV");
    CHECK(fmt::format("{:.3f}", in(1.0_mb, nb)) == "1000000.000 nb");
}
