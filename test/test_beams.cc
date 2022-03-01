#include "catch2/catch.hpp"

#include "nuchic/Beams.hh"
#include "yaml-cpp/yaml.h"

#include <iostream>

TEST_CASE("Spectrum Beam", "[Beams]") {
    SECTION("Parse headers") {
        SECTION("Parse Achilles header") {
            YAML::Node beam = YAML::Load("Histogram: flux/miniboone.dat");
            nuchic::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "Achilles");
            CHECK(spectrum.MinEnergy() == 0.0);
            CHECK(spectrum.MaxEnergy() == 4450.0);
        }
        SECTION("Parse MiniBooNE header") {
            YAML::Node beam = YAML::Load("Histogram: flux/miniboone_nu.dat");
            nuchic::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "MiniBooNE");
            CHECK(spectrum.MinEnergy() == 0.0);
            CHECK(spectrum.MaxEnergy() == 10000.0);
        }
        SECTION("Parse T2K header") {
            YAML::Node beam = YAML::Load("Histogram: flux/T2K_nu.dat");
            nuchic::Spectrum spectrum(beam);
            CHECK(spectrum.Format() == "T2K");
            CHECK(spectrum.MinEnergy() == 0.0);
            CHECK(spectrum.MaxEnergy() == 30000.0);
        }
    }

    SECTION("Parameters set correctly") {
        YAML::Node beam = YAML::Load("Histogram: flux/dummy.dat");
        nuchic::Spectrum spectrum(beam);

        CHECK(spectrum.NVariables() == 1);
        CHECK(spectrum.Flux({0.5}) == nuchic::FourVector(500, 0, 0, 500));
        std::vector<double> rans(1);
        CHECK(spectrum.GenerateWeight({500, 0, 0, 500}, rans) == 1./50.*1000);
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
        Energy: 200)beam"
        );
    
        auto beam = beams["Beams"].as<nuchic::Beam>();

        CHECK(beam.NBeams() == 2); 
        CHECK(beam.BeamIDs() == std::set<nuchic::PID>{12, 13});
        CHECK(beam.at(nuchic::PID(12)) == beam[nuchic::PID(12)]);
        CHECK(beam.Flux(nuchic::PID(12), {}) == nuchic::FourVector(100, 0, 0, 100));
        CHECK(beam.Flux(nuchic::PID(13), {}) == nuchic::FourVector(200, 0, 0, 200));
        std::vector<double> rans;
        CHECK(beam.GenerateWeight(nuchic::PID(12), {}, rans) == 1.0);
        CHECK(beam.GenerateWeight(nuchic::PID(13), {}, rans) == 1.0);
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
        Energy: 100)beam"
        );
    
        CHECK_THROWS_WITH(beams["Beams"].as<nuchic::Beam>(), "Multiple beams exist for PID: 12");
    }
}
