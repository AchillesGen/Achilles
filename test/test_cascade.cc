#include "catch2/catch.hpp"
#include "mock_classes.hh"

#include "nuchic/Cascade.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Interactions.hh"
#include "nuchic/Event.hh"

TEST_CASE("Initialize Cascade", "[Cascade]") {
    nuchic::Particles particles = {{nuchic::PID::proton()}, {nuchic::PID::neutron()}};

    SECTION("Kick Nucleon") {
        auto interaction = std::make_unique<MockInteraction>();
        auto nucleus = std::make_shared<MockNucleus>(); 

        REQUIRE_CALL(*nucleus, Nucleons())
            .LR_RETURN((particles))
            .TIMES(AT_LEAST(3));

        nuchic::Cascade cascade(std::move(interaction), nuchic::Cascade::ProbabilityType::Gaussian, nuchic::Cascade::InMedium::None);
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
        auto interaction = std::make_unique<MockInteraction>();
        auto nucleus = std::make_shared<MockNucleus>();
        std::shared_ptr<nuchic::Nucleus> tmp = nucleus;

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

        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
        cascade.Evolve(&event);

        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
        CHECK(hadrons[0].Radius() > radius);
    }

    SECTION("Evolve Nucleus") {
        auto interaction = std::make_unique<MockInteraction>();
        auto nucleus = std::make_shared<MockNucleus>();

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));

        REQUIRE_CALL(*nucleus, Radius())
            .TIMES(AT_LEAST(1))
            .RETURN(radius);

        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
        cascade.SetKicked(0);
        cascade.Evolve(nucleus);

        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
        CHECK(hadrons[0].Radius() > radius);
    }

    SECTION("NuWro Evolve") {
        auto interaction = std::make_unique<MockInteraction>();
        auto nucleus = std::make_shared<MockNucleus>();

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));

        REQUIRE_CALL(*nucleus, Radius())
            .TIMES(AT_LEAST(1))
            .RETURN(radius);

        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
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

    auto interaction = std::make_unique<MockInteraction>();
    auto nucleus = std::make_shared<MockNucleus>();
    spdlog::info("Mode: {}", mode);

    SECTION("Must have exactly one kicked") {
        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));

        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
        CHECK_THROWS_WITH(cascade.MeanFreePath(nucleus), "MeanFreePath: only one particle should be kicked.");

        cascade.SetKicked(0);
        cascade.SetKicked(1);
        CHECK_THROWS_WITH(cascade.MeanFreePath(nucleus), "MeanFreePath: only one particle should be kicked.");
    }

    SECTION("Must have internal test particle") {
        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(1)
            .LR_RETURN((hadrons));

        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
        cascade.SetKicked(1);
        CHECK_THROWS_WITH(cascade.MeanFreePath(nucleus),
            "MeanFreePath: kickNuc must have status -3 in order to accumulate DistanceTraveled.");
    }

    SECTION("Particle escapes marked correctly") {
        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));
        REQUIRE_CALL(*nucleus, Radius())
            .TIMES(AT_LEAST(1))
            .RETURN(radius);

        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
        cascade.SetKicked(0);
        CHECK_NOTHROW(cascade.MeanFreePath(nucleus));
        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
        CHECK(hadrons[0].Radius() > radius);
    }
}

TEST_CASE("NuWro Mean Free Path Mode", "[Cascade]") {

}
