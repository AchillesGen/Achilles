#include "catch2/catch.hpp"

#include "Achilles/InteractionHandler.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"

using achilles::PID;

TEST_CASE("Map ordering indpendent of key order", "[InteractionHandler]") {
    achilles::pid_compare compare;
    REQUIRE(compare({1, 2}, {2, 1}) == false);
}

TEST_CASE("Registering an interaction", "[InteractionHandler]") {
    achilles::InteractionHandler handler;
    achilles::ConstantInteraction interaction;
    interaction.AddInteraction({PID::proton(), PID::proton()},
                               {{{PID::proton(), PID::proton()}, 10}});
    interaction.AddInteraction({PID::proton(), PID::neutron()},
                               {{{PID::proton(), PID::neutron()}, 10}});
    interaction.AddInteraction({PID::neutron(), PID::neutron()},
                               {{{PID::neutron(), PID::neutron()}, 10}});
    handler.RegisterInteraction(interaction);
    REQUIRE(handler.RegisteredInteractions().size() == 3);
    REQUIRE(handler.EnabledChannel(PID::proton(), PID::proton()) == true);
    REQUIRE(handler.EnabledChannel(PID::proton(), PID::neutron()) == true);
    REQUIRE(handler.EnabledChannel(PID::neutron(), PID::proton()) == true);
    REQUIRE(handler.EnabledChannel(PID::neutron(), PID::neutron()) == true);
}

TEST_CASE("Testing Cross Section", "[InteractionHandler]") {
    achilles::InteractionHandler handler;
    achilles::ConstantInteraction interaction;
    interaction.AddInteraction({PID::proton(), PID::proton()},
                               {{{PID::proton(), PID::proton()}, 10}});

    handler.RegisterInteraction(interaction);
    achilles::Particle p1{PID::proton(), {}};
    achilles::Particle p2{PID::proton(), {}};
    auto results = handler.CrossSection(p1, p2);
    REQUIRE(results.size() == 1);
    REQUIRE(results[0].cross_section == 10);
    REQUIRE(results[0].particles.size() == 2);
    REQUIRE(results[0].particles[0] == PID::proton());
    REQUIRE(results[0].particles[1] == PID::proton());
}

TEST_CASE("Validate Interactions", "[InteractionHandler]") {
    achilles::InteractionHandler handler;
    achilles::ConstantInteraction interaction1, interaction2;
    interaction1.AddInteraction({PID::proton(), PID::proton()},
                                {{{PID::proton(), PID::proton()}, 10},
                                 {{PID::proton(), PID::neutron(), PID::pionp()}, 20}});
    interaction2.AddInteraction({PID::proton(), PID::proton()},
                                {{{PID::proton(), PID::proton()}, 10},
                                 {{PID::proton(), PID::neutron(), PID::pionp()}, 20}});

    handler.RegisterInteraction(interaction1);
    auto msg = fmt::format("Interaction already registered for particles: [2212, 2212]");
    REQUIRE_THROWS_MATCHES(handler.RegisterInteraction(interaction2), std::runtime_error,
                           Catch::Matchers::Message(msg));

    msg = fmt::format("Cross section not registered for particles: [211, 2212]");
    REQUIRE_THROWS_MATCHES(handler.CrossSection({PID::pionp(), {}}, {PID::proton(), {}}),
                           std::runtime_error, Catch::Matchers::Message(msg));

    msg = fmt::format("Phase space not registered for particles: [211, 2212]");
    REQUIRE_THROWS_MATCHES(
        handler.GenerateMomentum({PID::pionp(), {}}, {PID::proton(), {}}, {}, {}),
        std::runtime_error, Catch::Matchers::Message(msg));
}

TEST_CASE("Testing Multiple Final States", "[InteractionHandler]") {
    achilles::InteractionHandler handler;
    achilles::ConstantInteraction interaction;
    interaction.AddInteraction({PID::proton(), PID::proton()},
                               {{{PID::proton(), PID::proton()}, 10},
                                {{PID::proton(), PID::neutron(), PID::pionp()}, 20}});

    handler.RegisterInteraction(interaction);
    achilles::Particle p1{PID::proton(), {}};
    achilles::Particle p2{PID::proton(), {}};
    auto results = handler.CrossSection(p1, p2);
    REQUIRE(results.size() == 2);

    REQUIRE(results[0].cross_section == 10);
    REQUIRE(results[0].particles.size() == 2);
    REQUIRE(results[0].particles[0] == PID::proton());
    REQUIRE(results[0].particles[1] == PID::proton());

    REQUIRE(results[1].cross_section == 20);
    REQUIRE(results[1].particles.size() == 3);
    REQUIRE(results[1].particles[0] == PID::proton());
    REQUIRE(results[1].particles[1] == PID::neutron());
    REQUIRE(results[1].particles[2] == PID::pionp());
}

TEST_CASE("Channel Selection", "[InteractionHandler]") {
    achilles::InteractionHandler handler;
    achilles::ConstantInteraction interaction;
    interaction.AddInteraction({PID::proton(), PID::proton()},
                               {{{PID::proton(), PID::proton()}, 10},
                                {{PID::proton(), PID::neutron(), PID::pionp()}, 20}});

    handler.RegisterInteraction(interaction);
    achilles::Particle p1{PID::proton(), {}};
    achilles::Particle p2{PID::proton(), {}};
    auto results = handler.CrossSection(p1, p2);

    auto channel = handler.SelectChannel(results, 0.5);
    REQUIRE(channel.cross_section == 20);
    REQUIRE(channel.particles.size() == 3);
    REQUIRE(channel.particles[0] == PID::proton());
    REQUIRE(channel.particles[1] == PID::neutron());
    REQUIRE(channel.particles[2] == PID::pionp());

    REQUIRE_THROWS_MATCHES(handler.SelectChannel({}, 0.1), std::runtime_error,
                           Catch::Matchers::Message("No channels to select from"));
}
