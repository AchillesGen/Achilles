#include "catch2/catch.hpp"

#include "nuchic/Beams.hh"
#include "yaml-cpp/yaml.h"

#include <iostream>

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
        CHECK(beam.Flux(nuchic::PID(12), {}) == nuchic::FourVector(0, 0, 100, 100));
        CHECK(beam.Flux(nuchic::PID(13), {}) == nuchic::FourVector(0, 0, 200, 200));
        CHECK(beam.Weight(nuchic::PID(12), {0}) == 1.0);
        CHECK(beam.Weight(nuchic::PID(13), {0}) == 1.0);
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

    SECTION("Spectrum Beams") {
        YAML::Node beams = YAML::Load(R"beam(

Beams:
  - Beam:
      PID: 12
      Beam Params:
        Type: Spectrum
        Filename: dummy.txt)beam"
        );

        auto beam = beams["Beams"].as<nuchic::Beam>();
        auto spectrum = beam.at(nuchic::PID(12));

        CHECK_THROWS_WITH(spectrum -> NVariables(), "Spectrum Fluxes are not implemented");
        spdlog::info("Here");
        CHECK_THROWS_WITH(spectrum -> Flux({}), "Spectrum Fluxes are not implemented");
        spdlog::info("Here");
        CHECK_THROWS_WITH(spectrum -> Weight({}), "Spectrum Fluxes are not implemented");
    }
}
