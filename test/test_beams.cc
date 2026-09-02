// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"

#include "Achilles/Beams.hh"
#include "Achilles/YAML/Beams.hh"
#include "yaml-cpp/yaml.h"

#include <iostream>

namespace units = achilles::units;
using namespace achilles::units::literals;

TEST_CASE("Spectrum Beam", "[Beams]") {
    SECTION("Parse headers") {
        SECTION("Parse Achilles header") {
            YAML::Node beam = YAML::Load("Histogram: flux/miniboone.dat");
            achilles::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "Achilles");
            CHECK(spectrum.MinEnergy() == 0_MeV);
            CHECK(spectrum.MaxEnergy() == 4450_MeV);
        }
        SECTION("Parse MINVERvA Flux") {
            YAML::Node beam = YAML::Load("Histogram: flux/minerva_numu_fhc.dat");
            achilles::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "Achilles");
            CHECK(spectrum.MinEnergy() == 0_MeV);
            CHECK(spectrum.MaxEnergy() == 100_GeV);
        }
        SECTION("Parse MiniBooNE header") {
            YAML::Node beam = YAML::Load("Histogram: flux/miniboone_nu.dat");
            achilles::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "MiniBooNE");
            CHECK(spectrum.MinEnergy() == 0_MeV);
            CHECK(spectrum.MaxEnergy() == 10_GeV);
        }
        SECTION("Parse T2K header") {
            YAML::Node beam = YAML::Load("Histogram: flux/T2K_nu.dat");
            achilles::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "T2K");
            CHECK(spectrum.MinEnergy() == 0_MeV);
            CHECK(spectrum.MaxEnergy() == 30_GeV);
        }
    }

    SECTION("Parameters set correctly") {
        YAML::Node beam = YAML::Load("Histogram: flux/dummy.dat");
        achilles::Spectrum spectrum(beam);

        CHECK(spectrum.NVariables() == 1);
        CHECK(spectrum.Flux({0.5}, 0_MeV) == achilles::FourVector::FromNative({500, 0, 0, 500}));
        std::vector<double> rans(1);
        CHECK(spectrum.GenerateWeight({500_MeV, 0_MeV, 0_MeV, 500_MeV}, rans, 0_MeV) == 1.0);
        CHECK(rans[0] == 0.5);
    }
}

TEST_CASE("From YAML", "[Beams]") {
    SECTION("Multiple Monochromatic Beams") {
        YAML::Node beams = YAML::Load(R"beam(

Beams:
  - Beam:
      PID: 12
      Beam Params:
        Type: Monochromatic
        Energy: { value: 100, unit: MeV }
  - Beam:
      PID: 13
      Beam Params:
        Type: Monochromatic
        Energy: { value: 100, unit: MeV })beam");

        auto beam = beams["Beams"].as<achilles::Beam>();

        CHECK(beam.NBeams() == 2);
        CHECK(beam.BeamIDs() == std::set<achilles::PID>{12, 13});
        CHECK(beam.at(achilles::PID(12)) == beam[achilles::PID(12)]);
        CHECK(beam.Flux(achilles::PID(12), {}, 0_MeV) ==
              achilles::FourVector::FromNative({100, 0, 0, 100}));
        CHECK(beam.Flux(achilles::PID(13), {}, 0_MeV) ==
              achilles::FourVector::FromNative({100, 0, 0, 100}));
        std::vector<double> rans;
        CHECK(beam.GenerateWeight(achilles::PID(12), {}, rans, 0_MeV) == 1.0);
        CHECK(beam.GenerateWeight(achilles::PID(13), {}, rans, 0_MeV) == 1.0);
    }

    SECTION("Multiple Monochromatic Beams, throw on energy mismatch") {
        YAML::Node beams = YAML::Load(R"beam(

Beams:
  - Beam:
      PID: 12
      Beam Params:
        Type: Monochromatic
        Energy: { value: 100, unit: MeV }
  - Beam:
      PID: 13
      Beam Params:
        Type: Monochromatic
        Energy: { value: 200, unit: MeV })beam");

        CHECK_THROWS_WITH(beams["Beams"].as<achilles::Beam>(),
                          "Beams must have the same minimum energies. Got 100 MeV and 200 MeV");
    }

    SECTION("Throw on Identical Beams") {
        YAML::Node beams = YAML::Load(R"beam(

Beams:
  - Beam:
      PID: 12
      Beam Params:
        Type: Monochromatic
        Energy: { value: 100, unit: MeV }
  - Beam:
      PID: 12
      Beam Params:
        Type: Monochromatic
        Energy: { value: 100, unit: MeV })beam");

        CHECK_THROWS_WITH(beams["Beams"].as<achilles::Beam>(), "Multiple beams exist for PID: 12");
    }
}

TEST_CASE("FlatFlux", "[Beams]") {
    SECTION("Appropriately reads min and max energy") {
        YAML::Node flatflux = YAML::Load(R"beam(
MinEnergy: { value: 0, unit: MeV }
MaxEnergy: { value: 4000, unit: MeV }
)beam");

        achilles::FlatFlux flux(flatflux);
        CHECK(flux.MinEnergy() == 0_MeV);
        CHECK(flux.MaxEnergy() == 4000_MeV);
    }

    SECTION("Appropriately reads range") {
        YAML::Node flatflux = YAML::Load(R"beam(
Range: [{ value: 0, unit: MeV }, { value: 4000, unit: MeV }]
)beam");

        achilles::FlatFlux flux(flatflux);
        CHECK(flux.MinEnergy() == 0_MeV);
        CHECK(flux.MaxEnergy() == 4000_MeV);
    }

    SECTION("Throws if missing min/max energy or range") {
        YAML::Node flatflux = YAML::Load(R"beam(
)beam");

        CHECK_THROWS_AS(achilles::FlatFlux(flatflux), std::runtime_error);
    }

    SECTION("Throws if both min/max energy and range are given") {
        YAML::Node flatflux = YAML::Load(R"beam(
MinEnergy: { value: 0, unit: MeV }
MaxEnergy: { value: 4000, unit: MeV }
Range: [{ value: 0, unit: MeV }, { value: 4000, unit: MeV }]
)beam");

        CHECK_THROWS_AS(achilles::FlatFlux(flatflux), std::runtime_error);
    }
}
