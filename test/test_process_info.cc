#include "catch2/catch.hpp"

#include "Achilles/ProcessInfo.hh"
#include <sstream>

TEST_CASE("ProcessInfo", "[ProcessInfo]") {
    achilles::Process_Info info({achilles::PID::electron(), achilles::PID::electron()});
    info.state = {{achilles::PID::proton()}, {achilles::PID::proton()}};

    SECTION("Multiplicity is correct") {
        CHECK(info.Multiplicity() == 4);
    }

    SECTION("Masses are correct") {
        std::vector<double> expected{pow(achilles::ParticleInfo(achilles::PID::proton()).Mass(), 2),
                                     pow(achilles::ParticleInfo(achilles::PID::electron()).Mass(), 2)};

        CHECK(info.Masses() == expected);
    }

    SECTION("Ids are correct") {
        std::vector<long> expected{2212, 11, 2212, 11};

        CHECK(info.Ids() == expected);
    }

    SECTION("Output correct") {
        std::string expected = "Process_Info(PIDs = [11, 11], State = [[2212] -> [2212]])";
        std::stringstream output;
        output << info;

        CHECK(output.str() == expected);
    }
}

TEST_CASE("ProcessInfo YAML", "[ProcessInfo]") {
    YAML::Node node = YAML::Load(R"node(
        Model: SM
        Final States: [14, 11, -11]
    )node");

    auto info = node.as<achilles::Process_Info>();

    CHECK(info.ids == std::vector<achilles::PID>{14, 11, -11});
}
