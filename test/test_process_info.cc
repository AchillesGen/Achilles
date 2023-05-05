#include "catch2/catch.hpp"

#include "Achilles/ProcessInfo.hh"
#include <sstream>

TEST_CASE("ProcessInfo", "[ProcessInfo]") {
    achilles::ProcessInfo info({achilles::PID::electron(), {achilles::PID::electron()}});
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};

    SECTION("Multiplicity is correct") {
        CHECK(info.Multiplicity() == 4);
    }

    SECTION("Masses are correct") {
        std::vector<double> expected{
            pow(achilles::ParticleInfo(achilles::PID::electron()).Mass(), 2),
            pow(achilles::ParticleInfo(achilles::PID::proton()).Mass(), 2)};

        CHECK(info.Masses() == expected);
    }

    SECTION("Ids are correct") {
        std::vector<long> expected{11, 2212, 11, 2212};

        CHECK(info.Ids() == expected);
    }

    SECTION("Output correct") {
        std::string expected = "Process_Info([11, 2212] -> [11, 2212])";
        std::stringstream output;
        output << info;

        CHECK(output.str() == expected);
    }
}

TEST_CASE("ProcessInfo YAML", "[ProcessInfo]") {
    YAML::Node node = YAML::Load(R"node(
        Leptons: [14, [14, 11, -11]]
    )node");

    auto info = node.as<achilles::ProcessInfo>();

    CHECK(info.m_leptonic.first == achilles::PID{14});
    CHECK(info.m_leptonic.second == std::vector<achilles::PID>{14, 11, -11});
}
