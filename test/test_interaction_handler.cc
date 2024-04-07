#include "catch2/catch.hpp"

#include "Achilles/InteractionHandler.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"

TEST_CASE("Map ordering indpendent of key order", "[InteractionHandler]") {
    achilles::pid_compare compare;
    REQUIRE(compare({1, 2}, {2, 1}) == false);
}

TEST_CASE("Registering an interaction", "[InteractionHandler]") {
    achilles::InteractionHandler handler;
    achilles::ConstantInteraction interaction(10);
    handler.RegisterInteraction(interaction);
    REQUIRE(handler.RegisteredInteractions().size() == 3);
    REQUIRE(handler.EnabledChannel(achilles::PID::proton(), achilles::PID::proton()) == true);
    REQUIRE(handler.EnabledChannel(achilles::PID::proton(), achilles::PID::neutron()) == true);
    REQUIRE(handler.EnabledChannel(achilles::PID::neutron(), achilles::PID::proton()) == true);
    REQUIRE(handler.EnabledChannel(achilles::PID::neutron(), achilles::PID::neutron()) == true);
}

TEST_CASE("Testing Cross Section", "[InteractionHandler]") {
    achilles::InteractionHandler handler;
    achilles::ConstantInteraction interaction(10);
    handler.RegisterInteraction(interaction);
    achilles::Particle p1{achilles::PID::proton(), {}};
    achilles::Particle p2{achilles::PID::proton(), {}};
    auto results = handler.CrossSection(p1, p2);
    REQUIRE(results.size() == 1);
    REQUIRE(results[0].cross_section == 10);
    REQUIRE(results[0].particles.size() == 2);
    REQUIRE(results[0].particles[0] == achilles::PID::proton());
    REQUIRE(results[0].particles[1] == achilles::PID::proton());
}
