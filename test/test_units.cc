// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Tests for the strong-typed natural-units layer (Achilles/Units.hh).

#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_string.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Units.hh"
#include "Achilles/UnitsFormat.hh"
#include "Achilles/UnitsIO.hh"
#include "Achilles/UnitsSchema.hh"

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

// ---------------------------------------------------------------------------
// File-boundary enforcement: a bare number becomes a Quantity in exactly one
// place, and only by naming its unit.
// ---------------------------------------------------------------------------

TEST_CASE("Unit names are looked up per dimension", "[Units][IO]") {
    namespace io = achilles::units::io;

    CHECK(io::unit_from_string<1>("GeV").scale() == GeV.scale());
    CHECK(io::unit_from_string<-1>("fm").scale() == fm.scale());
    CHECK(io::unit_from_string<-2>("nb").scale() == nb.scale());

    SECTION("A unit of the wrong dimension is simply not in the table") {
        CHECK_THROWS_WITH(io::unit_from_string<1>("fm"),
                          Catch::Matchers::ContainsSubstring("not a known ENERGY unit"));
        CHECK_THROWS_WITH(io::unit_from_string<-1>("GeV"),
                          Catch::Matchers::ContainsSubstring("not a known LENGTH unit"));
        CHECK_THROWS_WITH(io::unit_from_string<-2>("MeV"),
                          Catch::Matchers::ContainsSubstring("not a known CROSS-SECTION unit"));
    }

    SECTION("An unknown unit throws rather than defaulting") {
        CHECK_THROWS_WITH(io::unit_from_string<1>("furlong"),
                          Catch::Matchers::ContainsSubstring("not a known ENERGY unit"));
        CHECK_THROWS_WITH(io::from_declared<1>(1.0, ""),
                          Catch::Matchers::ContainsSubstring("refusing to guess"));
    }
}

TEST_CASE("Declared, specified and range-checked values", "[Units][IO]") {
    namespace io = achilles::units::io;

    // (A) the file declares its unit
    CHECK(io::from_declared<1>(2.2, "GeV") == 2200_MeV);
    CHECK(io::from_declared<-1>(1.5, "fm") == 1.5_fm);

    // (B) the format fixes the unit -- LHE energies are GeV by the accord
    CHECK(io::from_lhe_energy(2.2) == 2.2_GeV);
    STATIC_REQUIRE(std::is_same<decltype(io::from_lhe_energy(1.0)), Energy>::value);

    // (C) range backstop
    CHECK(io::expect_range<1>(500_MeV, 0_MeV, 1000_MeV, "test") == 500_MeV);
    CHECK_THROWS_WITH(io::expect_range<1>(500_GeV, 0_MeV, 1000_MeV, "test"),
                      Catch::Matchers::ContainsSubstring("outside its physical range"));
}

TEST_CASE("Config quantities must name their unit", "[Units][IO]") {
    SECTION("A bare scalar is refused") {
        auto node = YAML::Load("beam: 30");
        CHECK_THROWS_WITH(node["beam"].as<Energy>(),
                          Catch::Matchers::ContainsSubstring("is a bare number"));
    }

    SECTION("A declared unit is applied") {
        auto node = YAML::Load("beam: { value: 30, unit: GeV }");
        CHECK(node["beam"].as<Energy>() == 30_GeV);
        CHECK(node["beam"].as<Energy>().in(MeV) == 30000.0);
    }

    SECTION("A wrong-dimension unit throws at load") {
        auto node = YAML::Load("beam: { value: 30, unit: fm }");
        CHECK_THROWS_WITH(node["beam"].as<Energy>(),
                          Catch::Matchers::ContainsSubstring("not a known ENERGY unit"));
        auto step = YAML::Load("step: { value: 0.04, unit: MeV }");
        CHECK_THROWS_WITH(step["step"].as<Length>(),
                          Catch::Matchers::ContainsSubstring("not a known LENGTH unit"));
    }

    SECTION("A half-written quantity throws") {
        auto no_unit = YAML::Load("beam: { value: 30 }");
        CHECK_THROWS(no_unit["beam"].as<Energy>());
        auto no_value = YAML::Load("beam: { unit: GeV }");
        CHECK_THROWS(no_value["beam"].as<Energy>());
    }

    SECTION("Round-trips through canonical units") {
        YAML::Node out;
        out["beam"] = 2.2_GeV;
        CHECK(out["beam"]["unit"].as<std::string>() == "MeV");
        CHECK(out["beam"].as<Energy>() == 2.2_GeV);
    }
}

TEST_CASE("Particle table columns declare their units once", "[Units][IO]") {
    namespace io = achilles::units::io;

    SECTION("An absent block means the documented MeV convention") {
        auto u = io::resolve_particle_units(YAML::Node());
        CHECK(io::particle_mass(938.27, u) == 938.27_MeV);
        CHECK(io::particle_width(0.0, u) == 0_MeV);
    }

    SECTION("A declared block overrides it and is dimension-checked") {
        auto node = YAML::Load("units: { mass: GeV, width: MeV }");
        auto u = io::resolve_particle_units(node["units"]);
        CHECK(io::particle_mass(0.93827, u).in(MeV) == Approx(938.27));
        CHECK(io::particle_width(1.0, u) == 1_MeV);

        auto bad = YAML::Load("units: { mass: fm }");
        CHECK_THROWS_WITH(io::resolve_particle_units(bad["units"]),
                          Catch::Matchers::ContainsSubstring("not a known ENERGY unit"));
    }

    SECTION("Range backstop catches an implausible mass") {
        auto node = YAML::Load("units: { mass: TeV }");
        auto u = io::resolve_particle_units(node["units"]);
        CHECK_THROWS_WITH(io::particle_mass(938.27, u),
                          Catch::Matchers::ContainsSubstring("outside its physical range"));
    }
}

TEST_CASE("The shipped Particles.yml loads in MeV", "[Units][IO]") {
    const achilles::ParticleInfo proton(achilles::PID::proton());
    const achilles::ParticleInfo muon(achilles::PID::muon());

    STATIC_REQUIRE(std::is_same<decltype(proton.Mass()), Energy>::value);
    CHECK(proton.Mass().in(MeV) == Approx(938.27).epsilon(1e-9));
    CHECK(muon.Mass().in(MeV) == Approx(105.7).epsilon(1e-3));
    CHECK(proton.Mass().in(GeV) == Approx(0.93827).epsilon(1e-9));
    CHECK(proton.Width() == 0_MeV);
}
