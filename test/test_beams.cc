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

    SECTION("Spectrum Beams") {
        YAML::Node beam = YAML::Load("Histogram: dummy.txt");
        nuchic::Spectrum spectrum(beam);

        for(size_t i = 0; i < 4400; ++i) {
            auto enu = static_cast<double>(i);
            nuchic::FourVector mom{enu, 0, 0, enu};
            std::vector<double> rans(1);
            double flux = spectrum.GenerateWeight(mom, rans);
            std::cout << fmt::format("{},{}\n", enu, flux);
        }
    }
}
