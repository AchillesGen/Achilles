#include "catch2/catch.hpp" 
#include "mock_classes.hh"

#include <sstream>

#include "nuchic/EventWriter.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Version.hh"

TEST_CASE("Builtin", "[EventWriter]") {
    SECTION("Write Header") {
        std::stringstream ss;
        nuchic::NuchicWriter writer(&ss);

        writer.WriteHeader("dummy.txt");
        
        std::string expected = fmt::format(R"result(Nuchic Version: {0}
{1:-^40}

{1:-^40}

)result", NUCHIC_VERSION, "");

        CHECK(ss.str() == expected);
    }

    SECTION("Writing Event") {
        std::stringstream ss;
        nuchic::NuchicWriter writer(&ss);

        static constexpr nuchic::FourVector hadron0{65.4247, 26.8702, -30.5306, -10.9449};
        static constexpr nuchic::FourVector hadron1{1560.42, -78.4858, -204.738, 1226.89};
        nuchic::Particles particles = {
            {nuchic::PID::proton(), hadron0, {}, nuchic::ParticleStatus::initial_state},
            {nuchic::PID::proton(), hadron1, {}, nuchic::ParticleStatus::final_state}};
        nuchic::NuclearRemnant remnant(11, 5);
        
        MockEvent event;
        REQUIRE_CALL(event, Particles())
            .TIMES(1)
            .LR_RETURN((particles));
        REQUIRE_CALL(event, Remnant())
            .TIMES(1)
            .LR_RETURN((remnant));
        REQUIRE_CALL(event, Weight())
            .TIMES(1)
            .RETURN(1.0);

        writer.Write(event);

        std::string expected = fmt::format(R"result(Event: 1
  Particles:
  - {}
  - {}
  - {}
  Weight: {}
)result", particles[0], particles[1], remnant, 1.0);

        CHECK(ss.str() == expected);
    }
}
