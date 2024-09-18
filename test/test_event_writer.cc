#include "catch2/catch.hpp"
#include "mock_classes.hh"

#include <sstream>

#include "Achilles/EventWriter.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Version.hh"

using achilles::Particle;

TEST_CASE("Builtin", "[EventWriter]") {
    SECTION("Write Header") {
        std::stringstream ss;
        achilles::AchillesWriter writer(&ss);
        std::vector<achilles::ProcessGroup> group;

        writer.WriteHeader("dummy.txt", group);

        std::string expected = fmt::format(R"result(Achilles Version: {0}
{1:-^40}

{1:-^40}

)result",
                                           ACHILLES_VERSION, "");

        CHECK(ss.str() == expected);
    }

    SECTION("Writing Event") {
        std::stringstream ss;
        achilles::AchillesWriter writer(&ss);

        static constexpr achilles::FourVector hadron0{65.4247, 26.8702, -30.5306, -10.9449};
        static constexpr achilles::FourVector hadron1{1560.42, -78.4858, -204.738, 1226.89};
        achilles::Particles particles = {
            Particle{achilles::PID::proton(), hadron0, {}, achilles::ParticleStatus::initial_state},
            Particle{achilles::PID::proton(), hadron1, {}, achilles::ParticleStatus::final_state}};
        achilles::NuclearRemnant remnant(11, 5);

        const MockEvent event{};
        double wgt = 1.0;
        REQUIRE_CALL(event, Particles()).TIMES(1).LR_RETURN((particles));
        REQUIRE_CALL(event, Remnant()).TIMES(1).LR_RETURN((remnant));
        REQUIRE_CALL(event, Weight()).TIMES(1).LR_RETURN((wgt));

        writer.Write(event);

        std::string expected = fmt::format(R"result(Event: 1
  Particles:
  - {}
  - {}
  - {}
  Weight: {}
)result",
                                           particles[0], particles[1], remnant, 1.0);

        CHECK(ss.str() == expected);
    }
}
