#include "catch2/catch.hpp"

#include "nuchic/CombinedCuts.hh"
#include "nuchic/Cuts.hh"
#include "nuchic/OneParticleCuts.hh"
#include "nuchic/TwoParticleCuts.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Units.hh"

#include "catch_utils.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include <iostream>

nuchic::FourVector const& RandomMomentumGenerator::get() const {
    return current_momentum;
}
// This helper function provides a nicer UX when instantiating the generator
// Notice that it returns an instance of GeneratorWrapper<std::string>, which
// is a value-wrapper around std::unique_ptr<IGenerator<std::string>>.
Catch::Generators::GeneratorWrapper<nuchic::FourVector> randomMomentum(double max) {
    return Catch::Generators::GeneratorWrapper<nuchic::FourVector>(std::unique_ptr<Catch::Generators::IGenerator<nuchic::FourVector>>(new RandomMomentumGenerator(max)));
}

template<class T>
struct TestCutBase : public nuchic::CutBase<T> {
    TestCutBase(YAML::Node node) : nuchic::CutBase<T>(node) {}
    bool MakeCut(const T& val) const { return this->CheckCut(val); }
};

TEMPLATE_TEST_CASE("Base Cut", "[Cuts]", double, float, size_t, int) {
    SECTION("Check for valid syntax") {
        YAML::Node node;
        REQUIRE_THROWS_AS(TestCutBase<TestType>(node), std::runtime_error);
        node["min"] = 30;
        REQUIRE_NOTHROW(TestCutBase<TestType>(node));
        node["max"] = 100;
        REQUIRE_NOTHROW(TestCutBase<TestType>(node));
        node["range"] = std::vector<TestType>{30, 50};
        REQUIRE_THROWS_AS(TestCutBase<TestType>(node), std::runtime_error);
        YAML::Node node2;
        node2["range"] = std::vector<TestType>{30, 50};
        REQUIRE_NOTHROW(TestCutBase<TestType>(node2));
    }

    SECTION("Cuts work as expected") {
        YAML::Node node;
        node["min"] = 10;
        auto cut = TestCutBase<TestType>(node);
        auto val = static_cast<TestType>(GENERATE(take(30, random(0., 100.))));
        CHECK((val > 10) == cut.MakeCut(val)); 
    }
}

TEST_CASE("Single particle cuts", "[Cuts]") {
    auto mom = GENERATE(take(30, randomMomentum(1000)));

    SECTION("Energy Cut") {
        YAML::Node node;
        node["min"] = 10;
        nuchic::EnergyCut cut(node);
        CHECK(cut.MakeCut(mom) == (mom.E() > 10));
    }

    SECTION("Momentum Cut") {
        YAML::Node node;
        node["min"] = 10;
        nuchic::MomentumCut cut(node);
        CHECK(cut.MakeCut(mom) == (mom.P() > 10));
    }

    SECTION("AngleTheta Cut") {
        YAML::Node node;
        node["min"] = 10;
        nuchic::AngleThetaCut cut(node);
        CHECK(cut.MakeCut(mom) == (mom.Theta()*180/M_PI > 10));
    }

    SECTION("Transverse Momentum Cut") {
        YAML::Node node;
        node["min"] = 10;
        nuchic::TransverseMomentumCut cut(node);
        CHECK(cut.MakeCut(mom) == (mom.Pt() > 10));
    }

    SECTION("ETheta2 Cut") {
        YAML::Node node;
        node["min"] = 10;
        nuchic::ETheta2Cut cut(node);
        CHECK(cut.MakeCut(mom) == (mom.E()*pow(mom.Theta(), 2) > 10));
    }
}

TEST_CASE("Two particle cuts", "[Cuts]") {
    auto mom1 = GENERATE(take(30, randomMomentum(1000)));
    auto mom2 = GENERATE(take(30, randomMomentum(1000)));

    SECTION("DeltaTheta Cut") {
        YAML::Node node;
        node["min"] = 10;
        nuchic::DeltaThetaCut cut(node);
        CHECK(cut.MakeCut(mom1, mom2) == (std::acos(mom1.CosAngle(mom2))*180/M_PI > 10));
    }

    SECTION("InvariantMass Cut") {
        YAML::Node node;
        node["min"] = 10;
        nuchic::InvariantMassCut cut(node);
        CHECK(cut.MakeCut(mom1, mom2) == ((mom1+mom2).M() > 10));
    }
}

TEST_CASE("CutCollection", "[Cuts]") {
    nuchic::CutCollection cuts;
    auto mom1 = GENERATE(take(30, randomMomentum(1000)));
    auto mom2 = GENERATE(take(30, randomMomentum(1000)));
    auto part1 = nuchic::Particle(nuchic::PID::electron(), mom1);
    auto part2 = nuchic::Particle(nuchic::PID::muon(), mom2);
    auto part3 = nuchic::Particle(nuchic::PID::electron(), mom2);
    auto part4 = nuchic::Particle(-nuchic::PID::electron(), mom2);

    SECTION("Single PIDs work") {
        YAML::Node node;
        node["min"] = 10;
        auto cut = nuchic::CutFactory<nuchic::OneParticleCut>::InitializeCut("Energy", node);
        cuts.AddCut({nuchic::PID::electron()}, std::move(cut)); 
        CHECK(cuts.EvaluateCuts({part1}) == (mom1.E() > 10));
        CHECK(cuts.EvaluateCuts({part2}) == true);
        CHECK(cuts.EvaluateCuts({part1, part2}) == (mom1.E() > 10));

        auto cut2 = nuchic::CutFactory<nuchic::TwoParticleCut>::InitializeCut("DeltaTheta", node);
        cuts.AddCut({nuchic::PID::electron()}, std::move(cut2));
        bool pass_cuts = (mom1.E() > 10 && mom2.E() > 10 && std::acos(mom1.CosAngle(mom2))*180/M_PI > 10);
        CHECK(cuts.EvaluateCuts({part1, part3}) == pass_cuts);
    }

    SECTION("Sets of PIDs work") {
        YAML::Node node;
        node["min"] = 10;
        auto cut = nuchic::CutFactory<nuchic::OneParticleCut>::InitializeCut("Energy", node);
        cuts.AddCut({nuchic::PID::electron(), nuchic::PID::muon()}, std::move(cut)); 
        CHECK(cuts.EvaluateCuts({part1, part2}) == (mom1.E() > 10 && mom2.E() > 10));

        auto cut2 = nuchic::CutFactory<nuchic::TwoParticleCut>::InitializeCut("DeltaTheta", node);
        cuts.AddCut({nuchic::PID::electron(), nuchic::PID::muon()}, std::move(cut2));
        bool pass_cuts = (mom1.E() > 10 && mom2.E() > 10 && std::acos(mom1.CosAngle(mom2))*180/M_PI > 10);
        CHECK(cuts.EvaluateCuts({part1, part2}) == pass_cuts);
    }

    SECTION("YAML correctly builds CutCollection") {
        YAML::Node node = YAML::Load(R"node(
        cuts:
            - Type: Energy
              PIDs: [11, -11]
              min: 10
            - Type: InvariantMass
              PIDs: [11, -11]
              range: [80, 100]
            - Type: AngleTheta
              PIDs: 13
              max: 30
            - Type: DeltaTheta
              PIDs: [11, -11]
              range: [[10, 30], [60, 80]]
        )node");

        cuts = node["cuts"].as<nuchic::CutCollection>();
        bool pass_cuts = (mom1.E() > 10) && (mom2.E() > 10);
        pass_cuts &= (mom1+mom2).M() > 80 && (mom1+mom2).M() < 100;
        double dtheta = std::acos(mom1.CosAngle(mom2))*180/M_PI;
        pass_cuts &= (dtheta > 10 && dtheta < 30) || (dtheta > 60 && dtheta < 80);
        CHECK(cuts.EvaluateCuts({part1, part4}) == pass_cuts);
    }
}
