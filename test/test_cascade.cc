#include "catch2/catch.hpp"
#include "catch2/trompeloeil.hpp"

#include "nuchic/Cascade.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"
#include "nuchic/InteractionComponent.hh"
#include "nuchic/Event.hh"

class MockNucleus : public trompeloeil::mock_interface<nuchic::Nucleus> {
    static constexpr bool trompeloeil_movable_mock = true;
    MAKE_MOCK0(Nucleons, nuchic::Particles&(), noexcept override);
    IMPLEMENT_MOCK0(GenerateConfig);
    MAKE_CONST_MOCK0(Radius, const double&(), noexcept override);
};

class MockInteraction : public trompeloeil::mock_interface<nuchic::InteractionComponent> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK2(CrossSection);
    IMPLEMENT_CONST_MOCK2(CrossSections);
    IMPLEMENT_CONST_MOCK2(GenerateFinalState);
    IMPLEMENT_CONST_MOCK0(InteractionType);
};

class MockEvent : public trompeloeil::mock_interface<nuchic::Event> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(CurrentNucleus);
    IMPLEMENT_MOCK0(Hadrons);
};

TEST_CASE("Initialize Cascade", "[Cascade]") {
    nuchic::Particles particles = {{nuchic::PID::proton()}, {nuchic::PID::neutron()}};

    SECTION("Kick Nucleon") {
        nuchic::CascadeInteraction interaction;
        auto nucleus = std::make_shared<MockNucleus>(); 
        auto interaction_component = std::make_unique<MockInteraction>();

        REQUIRE_CALL(*interaction_component, InteractionType())
            .LR_RETURN((InteractionComponentType::NucleonNucleon))
            .TIMES(1);

        REQUIRE_CALL(*nucleus, Nucleons())
            .LR_RETURN((particles))
            .TIMES(AT_LEAST(3));

        interaction.AddComponent(std::move(interaction_component));
        nuchic::Cascade cascade(std::move(interaction), nuchic::Cascade::ProbabilityType::Gaussian);
        cascade.Kick(nucleus, {100, 0, 0, 0}, {10, 0});
        CHECK(particles[0].Status() == nuchic::ParticleStatus::propagating);
        CHECK(particles[1].Status() == nuchic::ParticleStatus::background);
        particles[0].Status() = nuchic::ParticleStatus::background;

        cascade.Reset();
        cascade.Kick(nucleus, {100, 0, 0, 0}, {0, 10});
        CHECK(particles[0].Status() == nuchic::ParticleStatus::background);
        CHECK(particles[1].Status() == nuchic::ParticleStatus::propagating);
    }
}

TEST_CASE("Evolve States", "[Cascade]") {
    nuchic::Particles hadrons = {{nuchic::PID::proton(), {100, 0, 0, 1000},
                                 {0, 0, 0}, nuchic::ParticleStatus::propagating}};
    constexpr double radius = 1;

    auto mode = GENERATE(nuchic::Cascade::ProbabilityType::Gaussian,
                         nuchic::Cascade::ProbabilityType::Pion,
                         nuchic::Cascade::ProbabilityType::Cylinder);

    SECTION("Evolve Event") {
        MockEvent event;
        nuchic::CascadeInteraction interaction;
        auto nucleus = std::make_shared<MockNucleus>();
        std::shared_ptr<nuchic::Nucleus> tmp = nucleus;
        auto interaction_component = std::make_unique<MockInteraction>();

        REQUIRE_CALL(*interaction_component, InteractionType())
            .LR_RETURN((InteractionComponentType::NucleonNucleon))
            .TIMES(1);

        REQUIRE_CALL(event, Hadrons())
            .TIMES(AT_LEAST(2))
            .LR_RETURN((hadrons));
        REQUIRE_CALL(event, CurrentNucleus())
            .TIMES(1)
            .LR_RETURN((tmp));

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));

        REQUIRE_CALL(*nucleus, Radius())
            .TIMES(AT_LEAST(1))
            .RETURN(radius);

        interaction.AddComponent(std::move(interaction_component));
        nuchic::Cascade cascade(std::move(interaction), mode);
        cascade.Evolve(&event);

        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
        CHECK(hadrons[0].Radius() > radius);
    }

    SECTION("Evolve Nucleus") {
        nuchic::CascadeInteraction interaction;
        auto nucleus = std::make_shared<MockNucleus>();
        auto interaction_component = std::make_unique<MockInteraction>();

        REQUIRE_CALL(*interaction_component, InteractionType())
            .LR_RETURN((InteractionComponentType::NucleonNucleon))
            .TIMES(1);

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));

        REQUIRE_CALL(*nucleus, Radius())
            .TIMES(AT_LEAST(1))
            .RETURN(radius);

        interaction.AddComponent(std::move(interaction_component));
        nuchic::Cascade cascade(std::move(interaction), mode);
        cascade.SetKicked(0);
        cascade.Evolve(nucleus);

        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
        CHECK(hadrons[0].Radius() > radius);
    }

    SECTION("NuWro Evolve") {
        nuchic::CascadeInteraction interaction;
        auto nucleus = std::make_shared<MockNucleus>();
        auto interaction_component = std::make_unique<MockInteraction>();

        REQUIRE_CALL(*interaction_component, InteractionType())
            .LR_RETURN((InteractionComponentType::NucleonNucleon))
            .TIMES(1);

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));

        REQUIRE_CALL(*nucleus, Radius())
            .TIMES(AT_LEAST(1))
            .RETURN(radius);

        interaction.AddComponent(std::move(interaction_component));
        nuchic::Cascade cascade(std::move(interaction), mode);
        cascade.SetKicked(0);
        cascade.NuWro(nucleus);

        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
        CHECK(hadrons[0].Radius() > radius);
    }
}

TEST_CASE("Mean Free Path", "[Cascade]") {
    nuchic::Particles hadrons = {{nuchic::PID::proton(), {100, 0, 0, 1000},
                                 {0, 0, 0}, nuchic::ParticleStatus::internal_test},
                                 {nuchic::PID::proton(), {100, 0, 0, 1000},
                                 {0, 0, -1}, nuchic::ParticleStatus::background},
                                };
    constexpr double radius = 2;

    auto mode = GENERATE(nuchic::Cascade::ProbabilityType::Gaussian,
                         nuchic::Cascade::ProbabilityType::Pion,
                         nuchic::Cascade::ProbabilityType::Cylinder);

    nuchic::CascadeInteraction interaction;
    auto nucleus = std::make_shared<MockNucleus>();
    spdlog::info("Mode: {}", mode);

    SECTION("Must have exactly one kicked") {
        auto interaction_component = std::make_unique<MockInteraction>();

        REQUIRE_CALL(*interaction_component, InteractionType())
            .LR_RETURN((InteractionComponentType::NucleonNucleon))
            .TIMES(1);

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));

        interaction.AddComponent(std::move(interaction_component));
        nuchic::Cascade cascade(std::move(interaction), mode);
        CHECK_THROWS_WITH(cascade.MeanFreePath(nucleus), "MeanFreePath: only one particle should be kicked.");

        cascade.SetKicked(0);
        cascade.SetKicked(1);
        CHECK_THROWS_WITH(cascade.MeanFreePath(nucleus), "MeanFreePath: only one particle should be kicked.");
    }

    SECTION("Must have internal test particle") {
        auto interaction_component = std::make_unique<MockInteraction>();

        REQUIRE_CALL(*interaction_component, InteractionType())
            .LR_RETURN((InteractionComponentType::NucleonNucleon))
            .TIMES(1);

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(1)
            .LR_RETURN((hadrons));

        interaction.AddComponent(std::move(interaction_component));
        nuchic::Cascade cascade(std::move(interaction), mode);
        cascade.SetKicked(1);
        CHECK_THROWS_WITH(cascade.MeanFreePath(nucleus),
            "MeanFreePath: kickNuc must have status -3 in order to accumulate DistanceTraveled.");
    }

    SECTION("Particle escapes marked correctly") {
        auto interaction_component = std::make_unique<MockInteraction>();

        REQUIRE_CALL(*interaction_component, InteractionType())
            .LR_RETURN((InteractionComponentType::NucleonNucleon))
            .TIMES(1);

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));
        REQUIRE_CALL(*nucleus, Radius())
            .TIMES(AT_LEAST(1))
            .RETURN(radius);

        interaction.AddComponent(std::move(interaction_component));
        nuchic::Cascade cascade(std::move(interaction), mode);
        cascade.SetKicked(0);
        CHECK_NOTHROW(cascade.MeanFreePath(nucleus));
        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
        CHECK(hadrons[0].Radius() > radius);
    }
}

TEST_CASE("NuWro Mean Free Path Mode", "[Cascade]") {

}
