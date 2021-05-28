#include "catch2/catch.hpp"

#include "nuchic/Cuts.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Particle.hh"

#include "catch_utils.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include <iostream>

TEST_CASE("Single particle cuts", "[Cuts]") {
    SECTION("Cuts must have a valid type") {
        auto type = GENERATE(values({
                nuchic::CutType::Q2,
                nuchic::CutType::Energy,
                nuchic::CutType::Momentum,
                nuchic::CutType::InvariantMass,
                nuchic::CutType::TransverseMomentum,
                nuchic::CutType::AngleTheta,
                nuchic::CutType::AnglePhi,
                nuchic::CutType::ETheta2
            }));

        nuchic::Cut cut;
        CHECK_NOTHROW(cut.Add(type, 10));
    }

    SECTION("Single values define minimum cut") {
        nuchic::Cut cut;
        cut.Add(nuchic::CutType::Energy, 100);

        auto mom = GENERATE(take(30, randomMomentum(1000)));
        CHECK(cut(mom) == (mom.E() > 100));
    }

    SECTION("Single range expects value to be within range") {
        nuchic::Cut cut;
        cut.Add(nuchic::CutType::Energy, {{50, 100}});

        auto mom = GENERATE(take(30, randomMomentum(1000)));
        CHECK(cut(mom) == (mom.E() > 50 && mom.E() < 100));
    }

    SECTION("Negative number in range means no upper cut") {
        nuchic::Cut cut;
        cut.Add(nuchic::CutType::Energy, {{100, -1}});

        auto mom = GENERATE(take(30, randomMomentum(1000)));
        CHECK(cut(mom) == mom.E() > 100);
    }

    SECTION("Multiple ranges define multiple allowed regions") {
        nuchic::Cut cut;
        cut.Add(nuchic::CutType::Energy, {{50, 100}, {150, std::numeric_limits<double>::infinity()}});

        auto mom = GENERATE(take(30, randomMomentum(1000)));
        CHECK(cut(mom) == ((mom.E() > 50 && mom.E() < 100)
                           || (mom.E() > 150)));
    }

    SECTION("Can only have one cut of each allowed type except custom") {
        nuchic::Cut cut;
        cut.Add(nuchic::CutType::Energy, {{50, 100}});
        CHECK_THROWS_WITH(cut.Add(nuchic::CutType::Energy, -1),
                          "Cut: May only have one cut of each type except Custom. Found additional cut of type: nuchic::CutType::Energy");

        cut.Add([](const nuchic::FourVector&) { return true; }, -1);
        CHECK_NOTHROW(cut.Add([](const nuchic::FourVector&) { return true; }, -1));
    }
}

TEST_CASE("Single particle cuts from YAML", "[Cuts]") {
    YAML::Node node = YAML::Load(R"node(
    cuts:
      - cut:
          type: AngleTheta
      - cut:
          type: Energy
          range: [50, 100]
      - cut:
          type: Energy
          min: 50
          max: 100
      - cut:
          type: Energy
          max: 100
      - cut:
          type: Energy
          min: 10
      - cut:
          type: Energy
          range: [[50, 100], [200, 300]])node");

    SECTION("Cuts must have either a range or min / max option") {
        CHECK_THROWS(node["cuts"][0]["cut"].as<nuchic::Cut>());
    }

    SECTION("Cuts must have a valid type") {
        {
            auto cutNode = GENERATE(as<std::string>{},
                    "cut:\n    min: 10\n    type: Q2",
                    "cut:\n    min: 10\n    type: Energy",
                    "cut:\n    min: 10\n    type: Momentum",
                    "cut:\n    min: 10\n    type: Mass",
                    "cut:\n    min: 10\n    type: TransverseMomentum",
                    "cut:\n    min: 10\n    type: Theta",
                    "cut:\n    min: 10\n    type: Phi",
                    "cut:\n    min: 10\n    type: ETheta2"
                );
            auto yamlCut = YAML::Load(cutNode);

            CHECK_NOTHROW(yamlCut["cut"].as<nuchic::Cut>());
        }
    }

    SECTION("Cuts accept only minimum value") {
        CHECK_NOTHROW(node["cuts"][4]["cut"].as<nuchic::Cut>());

        auto cut = node["cuts"][4]["cut"].as<nuchic::Cut>();
        auto mom = GENERATE(take(30, randomMomentum(1000)));
        CHECK(cut(mom) == (mom.E() > 10));
    }

    SECTION("Cuts accept only a maximum value") {
        CHECK_NOTHROW(node["cuts"][3]["cut"].as<nuchic::Cut>());

        auto cut = node["cuts"][3]["cut"].as<nuchic::Cut>();
        auto mom = GENERATE(take(30, randomMomentum(1000)));
        CHECK(cut(mom) == (mom.E() < 100));
    }

    SECTION("Cuts accept both a minimum and maximum value") {
        CHECK_NOTHROW(node["cuts"][2]["cut"].as<nuchic::Cut>());

        auto cut = node["cuts"][2]["cut"].as<nuchic::Cut>();
        auto mom = GENERATE(take(30, randomMomentum(1000)));
        CHECK(cut(mom) == (mom.E() < 100 && mom.E() > 50));
    }

    SECTION("Cuts accept a range of allowed values") {
        CHECK_NOTHROW(node["cuts"][1]["cut"].as<nuchic::Cut>());

        auto cut = node["cuts"][1]["cut"].as<nuchic::Cut>();
        auto mom = GENERATE(take(30, randomMomentum(1000)));
        CHECK(cut(mom) == (mom.E() < 100 && mom.E() > 50));
    }

    SECTION("Cuts accept multiple ranges") {
        CHECK_NOTHROW(node["cuts"][5]["cut"].as<nuchic::Cut>());

        auto cut = node["cuts"][5]["cut"].as<nuchic::Cut>();
        auto mom = GENERATE(take(30, randomMomentum(1000)));
        CHECK(cut(mom) == ((mom.E() < 100 && mom.E() > 50)
                           || (mom.E() > 200 && mom.E() < 300)));

    }
}

TEST_CASE("Load all cuts into a map", "[Cuts]") {
    YAML::Node node = YAML::Load(R"node(
Cuts:
  e-:
      - cut:
          type: Energy
          range: [600, 900]
      - cut:
          type: Momentum
          min: 600
          max: 900
  12:
      - cut:
          type: Energy
          range: [200, 500]
      - cut:
          type: Momentum
          min: 200
          max: 500
    )node");

    SECTION("Load from a PID or a name") {
        CHECK_NOTHROW(node["Cuts"].as<nuchic::Cuts>());

        auto cuts = node["Cuts"].as<nuchic::Cuts>();

        auto mom = GENERATE(take(30, randomMomentum(1000)));

        // From name
        bool electronCuts = (mom.E() > 600 && mom.E() < 900);
        electronCuts &= (mom.P() > 600 && mom.P() < 900);
        CHECK(cuts[{nuchic::PID::electron()}](mom) == electronCuts);

        // From PID
        bool neutrinoCuts = (mom.E() > 200 && mom.E() < 500);
        neutrinoCuts &= (mom.P() > 200 && mom.P() < 500);
        CHECK(cuts[{nuchic::PID::nu_electron()}](mom) == neutrinoCuts);
    }
}
