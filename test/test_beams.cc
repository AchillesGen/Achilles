#include "catch2/catch.hpp"

#include "Achilles/Beams.hh"
#include "Achilles/Factory.hh"
#include "Achilles/FluxLoaders.hh"
#include "Achilles/YAML/Beams.hh"
#include "yaml-cpp/yaml.h"

#include <iostream>

TEST_CASE("Flux file loaders", "[Beams]") {
    SECTION("HistogramFluxLoader detects and parses the Achilles format") {
        YAML::Node node = YAML::Load("Histogram: flux/miniboone.dat");
        auto data = achilles::HistogramFluxLoader{}.Load(node);
        CHECK(data.format == "Achilles");
        CHECK(data.min_energy == 0.0);
        CHECK(data.max_energy == 4450.0);
        CHECK(data.bin_centers.size() == data.heights.size());
        CHECK(data.bin_centers.size() >= 2);
    }

    SECTION("HistogramFluxLoader detects and parses the MiniBooNE format") {
        YAML::Node node = YAML::Load("Histogram: flux/miniboone_nu.dat");
        auto data = achilles::HistogramFluxLoader{}.Load(node);
        CHECK(data.format == "MiniBooNE");
        CHECK(data.min_energy == 0.0);
        CHECK(data.max_energy == 10.0);
    }

    SECTION("HistogramFluxLoader detects and parses the T2K format") {
        YAML::Node node = YAML::Load("Histogram: flux/T2K_nu.dat");
        auto data = achilles::HistogramFluxLoader{}.Load(node);
        CHECK(data.format == "T2K");
        CHECK(data.min_energy == 0.0);
        CHECK(data.max_energy == 30.0);
    }
}

TEST_CASE("Flux factory dispatch", "[Beams]") {
    using FluxFactory = achilles::Factory<achilles::FluxType, const YAML::Node &>;

    SECTION("Known flux types self-register") {
        CHECK(FluxFactory::IsRegistered("Monochromatic"));
        CHECK(FluxFactory::IsRegistered("Spectrum"));
        CHECK(FluxFactory::IsRegistered("FlatFlux"));
        CHECK(FluxFactory::IsRegistered("PDFBeam"));
    }

    SECTION("Unknown flux type is not registered") {
        CHECK_FALSE(FluxFactory::IsRegistered("Nonexistent"));
    }

    SECTION("Decode dispatches through the factory") {
        YAML::Node node = YAML::Load("Type: FlatFlux\nMinEnergy: 100\nMaxEnergy: 500");
        auto flux = node.as<std::shared_ptr<achilles::FluxType>>();
        REQUIRE(flux != nullptr);
        CHECK(flux->Type() == "FlatFlux");
        CHECK(flux->MinEnergy() == 100);
        CHECK(flux->MaxEnergy() == 500);
    }

    SECTION("Decode of an unknown flux type throws") {
        YAML::Node node = YAML::Load("Type: Nonexistent");
        CHECK_THROWS(node.as<std::shared_ptr<achilles::FluxType>>());
    }
}

TEST_CASE("Spectrum Beam", "[Beams]") {
    SECTION("Parse headers") {
        SECTION("Parse Achilles header") {
            YAML::Node beam = YAML::Load("Histogram: flux/miniboone.dat");
            achilles::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "Achilles");
            CHECK(spectrum.MinEnergy() == 0.0);
            CHECK(spectrum.MaxEnergy() == 4450.0);
        }
        SECTION("Parse MINVERvA Flux") {
            YAML::Node beam = YAML::Load("Histogram: flux/minerva_numu_fhc.dat");
            achilles::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "Achilles");
            CHECK(spectrum.MinEnergy() == 0.0);
            CHECK(spectrum.MaxEnergy() == 100.0);
        }
        SECTION("Parse MiniBooNE header") {
            YAML::Node beam = YAML::Load("Histogram: flux/miniboone_nu.dat");
            achilles::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "MiniBooNE");
            CHECK(spectrum.MinEnergy() == 0.0);
            CHECK(spectrum.MaxEnergy() == 10.0);
        }
        SECTION("Parse T2K header") {
            YAML::Node beam = YAML::Load("Histogram: flux/T2K_nu.dat");
            achilles::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "T2K");
            CHECK(spectrum.MinEnergy() == 0.0);
            CHECK(spectrum.MaxEnergy() == 30.0);
        }
    }

    SECTION("Parameters set correctly") {
        YAML::Node beam = YAML::Load("Histogram: flux/dummy.dat");
        achilles::Spectrum spectrum(beam);

        CHECK(spectrum.NVariables() == 1);
        CHECK(spectrum.Flux({0.5}, 0) == achilles::FourVector(500, 0, 0, 500));
        std::vector<double> rans(1);
        CHECK(spectrum.GenerateWeight({500, 0, 0, 500}, rans, 0) == 1.0);
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
        Energy: 100
  - Beam:
      PID: 13
      Beam Params:
        Type: Monochromatic
        Energy: 100)beam");

        auto beam = beams["Beams"].as<achilles::Beam>();

        CHECK(beam.NBeams() == 2);
        CHECK(beam.BeamIDs() == std::set<achilles::PID>{12, 13});
        CHECK(beam.at(achilles::PID(12)) == beam[achilles::PID(12)]);
        CHECK(beam.Flux(achilles::PID(12), {}, 0) == achilles::FourVector(100, 0, 0, 100));
        CHECK(beam.Flux(achilles::PID(13), {}, 0) == achilles::FourVector(100, 0, 0, 100));
        std::vector<double> rans;
        CHECK(beam.GenerateWeight(achilles::PID(12), {}, rans, 0) == 1.0);
        CHECK(beam.GenerateWeight(achilles::PID(13), {}, rans, 0) == 1.0);
    }

    SECTION("Multiple Monochromatic Beams, throw on energy mismatch") {
        YAML::Node beams = YAML::Load(R"beam(

Beams:
  - Beam:
      PID: 12
      Beam Params:
        Type: Monochromatic
        Energy: 100
  - Beam:
      PID: 13
      Beam Params:
        Type: Monochromatic
        Energy: 200)beam");

        CHECK_THROWS_WITH(beams["Beams"].as<achilles::Beam>(),
                          "Beams must have the same minimum energies. Got 100 and 200");
    }

    SECTION("Throw on Identical Beams") {
        YAML::Node beams = YAML::Load(R"beam(

Beams:
  - Beam:
      PID: 12
      Beam Params:
        Type: Monochromatic
        Energy: 100
  - Beam:
      PID: 12
      Beam Params:
        Type: Monochromatic
        Energy: 100)beam");

        CHECK_THROWS_WITH(beams["Beams"].as<achilles::Beam>(), "Multiple beams exist for PID: 12");
    }
}
