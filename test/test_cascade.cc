#include "catch2/catch.hpp"
#include "mock_classes.hh"

#include "Achilles/Cascade.hh"
#include "Achilles/Event.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

using achilles::Particle;

TEST_CASE("Initialize Cascade", "[Cascade]") {
    achilles::Particles particles = {Particle{achilles::PID::proton()},
                                     Particle{achilles::PID::neutron()}};

    SECTION("Kick Nucleon") {
        achilles::InteractionHandler interaction;
        MockNucleus nucleus;
        MockEvent event;

        REQUIRE_CALL(event, Hadrons()).LR_RETURN((particles)).TIMES(AT_LEAST(3));

        achilles::Cascade cascade(
            std::move(interaction), achilles::Cascade::ProbabilityType::Gaussian,
            achilles::Cascade::Algorithm::Base, achilles::Cascade::InMedium::None);
        cascade.Kick(event, {0, 100, 0, 0}, {10, 0});
        CHECK(particles[0].Status() == achilles::ParticleStatus::propagating);
        CHECK(particles[1].Status() == achilles::ParticleStatus::background);
        particles[0].Status() = achilles::ParticleStatus::background;

        cascade.Reset();
        cascade.Kick(event, {0, 100, 0, 0}, {0, 10});
        CHECK(particles[0].Status() == achilles::ParticleStatus::background);
        CHECK(particles[1].Status() == achilles::ParticleStatus::propagating);
    }
}

TEST_CASE("Evolve States: 1 nucleon", "[Cascade]") {
    achilles::Particles hadrons = {Particle{achilles::PID::proton(),
                                            {1000, 100, 0, 0},
                                            {0, 0, 0},
                                            achilles::ParticleStatus::propagating}};
    constexpr double radius = 1;

    auto mode = GENERATE(achilles::Cascade::ProbabilityType::Gaussian,
                         achilles::Cascade::ProbabilityType::Pion,
                         achilles::Cascade::ProbabilityType::Cylinder);

    SECTION("Captured") {
        achilles::InteractionHandler interaction;
        MockNucleus nucleus;
        MockEvent event;
        auto potential = std::make_shared<MockPotential>();

        REQUIRE_CALL(event, Hadrons()).TIMES(AT_LEAST(2)).LR_RETURN((hadrons));
        REQUIRE_CALL(nucleus, GetPotential()).TIMES(AT_LEAST(1)).LR_RETURN((potential));

        REQUIRE_CALL(*potential, Hamiltonian(hadrons[0].Momentum().P(), hadrons[0].Position().P()))
            .TIMES(AT_LEAST(1))
            .RETURN(5);

        achilles::Cascade cascade(std::move(interaction), mode, achilles::Cascade::Algorithm::Base,
                                  achilles::Cascade::InMedium::None, true);
        cascade.SetKicked(0);
        cascade.Evolve(event, &nucleus);

        CHECK(hadrons[0].Status() == achilles::ParticleStatus::captured);
    }

    SECTION("Large Formation Zone") {
        achilles::InteractionHandler interaction;
        MockNucleus nucleus;
        MockEvent event;
        auto potential = std::make_shared<MockPotential>();

        REQUIRE_CALL(event, Hadrons()).TIMES(AT_LEAST(2)).LR_RETURN((hadrons));

        REQUIRE_CALL(nucleus, GetPotential()).TIMES(AT_LEAST(1)).RETURN(nullptr);
        REQUIRE_CALL(nucleus, Radius()).TIMES(AT_LEAST(1)).RETURN(radius);

        hadrons[0].SetFormationZone({10000, 0, 0, 0}, {88.2, 0, 0, 0});
        achilles::Cascade cascade(std::move(interaction), mode, achilles::Cascade::Algorithm::Base,
                                  achilles::Cascade::InMedium::None);
        cascade.SetKicked(0);
        cascade.Evolve(event, &nucleus);

        CHECK(hadrons[0].Status() == achilles::ParticleStatus::final_state);
    }

    SECTION("Evolve Event") {
        MockEvent event;
        achilles::InteractionHandler interaction;
        MockNucleus nucleus;

        REQUIRE_CALL(event, Hadrons()).TIMES(AT_LEAST(2)).LR_RETURN((hadrons));

        REQUIRE_CALL(nucleus, GetPotential()).TIMES(AT_LEAST(1)).RETURN(nullptr);
        REQUIRE_CALL(nucleus, Radius()).TIMES(AT_LEAST(1)).RETURN(radius);

        achilles::Cascade cascade(std::move(interaction), mode, achilles::Cascade::Algorithm::Base,
                                  achilles::Cascade::InMedium::None);
        cascade.Evolve(event, &nucleus);

        CHECK(hadrons[0].Status() == achilles::ParticleStatus::final_state);
        CHECK(hadrons[0].Radius() > radius);
    }

    // TODO: Restore MFP Approach
    // SECTION("NuWro Evolve") {
    //     achilles::InteractionHandler interaction;
    //     MockNucleus nucleus;
    //     MockEvent event;

    //     REQUIRE_CALL(event, Hadrons()).TIMES(AT_LEAST(2)).LR_RETURN((hadrons));

    //     REQUIRE_CALL(nucleus, GetPotential()).TIMES(AT_LEAST(1)).RETURN(nullptr);
    //     REQUIRE_CALL(nucleus, Radius()).TIMES(AT_LEAST(1)).RETURN(radius);

    //     achilles::Cascade cascade(std::move(interaction), mode,
    //     achilles::Cascade::Algorithm::MFP,
    //                               achilles::Cascade::InMedium::None);
    //     cascade.SetKicked(0);
    //     cascade.Evolve(event, &nucleus);

    //     CHECK(hadrons[0].Status() == achilles::ParticleStatus::final_state);
    //     CHECK(hadrons[0].Radius() > radius);
    // }
}

TEST_CASE("Evolve States: 3 nucleons", "[Cascade]") {
    achilles::Particles hadrons = {Particle{achilles::PID::proton(),
                                            {1000, 100, 0, 0},
                                            {0, 0, 0},
                                            achilles::ParticleStatus::propagating},
                                   Particle{achilles::PID::proton(),
                                            {achilles::Constant::mN, 0, 0, 0},
                                            {0.5, 0, 0},
                                            achilles::ParticleStatus::background},
                                   Particle{achilles::PID::neutron(),
                                            {achilles::Constant::mN, 0, 0, 0},
                                            {3, 0, 0},
                                            achilles::ParticleStatus::background}};
    constexpr double radius = 4;

    auto mode = GENERATE(achilles::Cascade::ProbabilityType::Gaussian,
                         achilles::Cascade::ProbabilityType::Pion,
                         achilles::Cascade::ProbabilityType::Cylinder);

    SECTION("Captured") {
        achilles::InteractionHandler interaction;
        MockEvent event;
        MockNucleus nucleus;
        auto potential = std::make_shared<MockPotential>();

        REQUIRE_CALL(event, Hadrons()).TIMES(AT_LEAST(2)).LR_RETURN((hadrons));

        REQUIRE_CALL(nucleus, GetPotential()).TIMES(AT_LEAST(1)).LR_RETURN((potential));
        REQUIRE_CALL(*potential, Hamiltonian(hadrons[0].Momentum().P(), hadrons[0].Position().P()))
            .TIMES(AT_LEAST(1))
            .RETURN(5);

        achilles::Cascade cascade(std::move(interaction), mode, achilles::Cascade::Algorithm::Base,
                                  achilles::Cascade::InMedium::None, true);
        cascade.SetKicked(0);
        cascade.Evolve(event, &nucleus);

        CHECK(hadrons[0].Status() == achilles::ParticleStatus::captured);
        CHECK(hadrons[1].Status() == achilles::ParticleStatus::background);
        CHECK(hadrons[2].Status() == achilles::ParticleStatus::background);
    }

    SECTION("Large Formation Zone") {
        achilles::InteractionHandler interaction;
        MockEvent event;
        MockNucleus nucleus;
        auto potential = std::make_shared<MockPotential>();

        REQUIRE_CALL(event, Hadrons()).TIMES(AT_LEAST(2)).LR_RETURN((hadrons));

        REQUIRE_CALL(nucleus, GetPotential()).TIMES(AT_LEAST(1)).RETURN(nullptr);
        REQUIRE_CALL(nucleus, Radius()).TIMES(AT_LEAST(1)).RETURN(radius);

        hadrons[0].SetFormationZone({10000, 0, 0, 0}, {88.2, 0, 0, 0});
        achilles::Cascade cascade(std::move(interaction), mode, achilles::Cascade::Algorithm::Base,
                                  achilles::Cascade::InMedium::None);
        cascade.SetKicked(0);
        cascade.Evolve(event, &nucleus);

        CHECK(hadrons[0].Status() == achilles::ParticleStatus::final_state);
        CHECK(hadrons[1].Status() == achilles::ParticleStatus::background);
        CHECK(hadrons[2].Status() == achilles::ParticleStatus::background);
    }

    // SECTION("Evolve Event") {
    //     MockEvent event;
    //     achilles::InteractionHandler interaction;
    //     auto nucleus = std::make_shared<MockNucleus>();
    //     std::shared_ptr<achilles::Nucleus> tmp = nucleus;

    //     REQUIRE_CALL(event, Hadrons()).TIMES(AT_LEAST(2)).LR_RETURN((hadrons));
    //     REQUIRE_CALL(event, CurrentNucleus()).TIMES(1).LR_RETURN((tmp));

    //     REQUIRE_CALL(*nucleus, Nucleons()).TIMES(AT_LEAST(2)).LR_RETURN((hadrons));
    //     REQUIRE_CALL(*nucleus, GetPotential()).TIMES(AT_LEAST(1)).RETURN(nullptr);
    //     REQUIRE_CALL(*nucleus, Rho(trompeloeil::gt(0))).TIMES(4).RETURN(0);

    //     REQUIRE_CALL(*nucleus, Radius()).TIMES(AT_LEAST(1)).RETURN(radius);

    //     std::pair<achilles::FourVector, achilles::FourVector> output{{1000, 80, 0, 0},
    //                                                                  {150, -20, 0, 0}};
    //     // REQUIRE_CALL(*interaction, CrossSection(trompeloeil::_,
    //     // hadrons[1])).TIMES(1).RETURN(1000); REQUIRE_CALL(*interaction,
    //     // CrossSection(trompeloeil::_, hadrons[2])).TIMES(1).RETURN(1000);
    //     // REQUIRE_CALL(*interaction, FinalizeMomentum(trompeloeil::_, trompeloeil::_,
    //     // trompeloeil::_))
    //     //     .TIMES(2)
    //     //     .LR_RETURN((output));

    //     achilles::Cascade cascade(interaction, mode, achilles::Cascade::Algorithm::Base,
    //                               achilles::Cascade::InMedium::None);
    //     cascade.Evolve(&event);

    //     CHECK(hadrons[0].Status() == achilles::ParticleStatus::final_state);
    //     CHECK(hadrons[1].Status() == achilles::ParticleStatus::background);
    //     CHECK(hadrons[2].Status() == achilles::ParticleStatus::background);
    //     CHECK(hadrons[0].Radius() > radius);
    // }

    // SECTION("NuWro Evolve") {
    //     auto interaction = std::make_unique<MockInteraction>();
    //     auto nucleus = std::make_shared<MockNucleus>();

    //     REQUIRE_CALL(*nucleus, Nucleons())
    //         .TIMES(2)
    //         .LR_RETURN((hadrons));

    //     REQUIRE_CALL(*nucleus, Radius())
    //         .TIMES(AT_LEAST(1))
    //         .RETURN(radius);

    //     achilles::Cascade cascade(std::move(interaction), mode,
    //     achilles::Cascade::InMedium::None); cascade.SetKicked(0); cascade.NuWro(nucleus);

    //     CHECK(hadrons[0].Status() == achilles::ParticleStatus::final_state);
    //     CHECK(hadrons[0].Radius() > radius);
    // }
}

TEST_CASE("Mean Free Path", "[Cascade]") {
    achilles::Particles hadrons = {
        Particle{achilles::PID::proton(),
                 {100, 0, 0, 1000},
                 {0, 0, 0},
                 achilles::ParticleStatus::internal_test},
        Particle{achilles::PID::proton(),
                 {100, 0, 0, 1000},
                 {0, 0, -1},
                 achilles::ParticleStatus::background},
    };
    constexpr double radius = 2;

    auto mode = GENERATE(achilles::Cascade::ProbabilityType::Gaussian,
                         achilles::Cascade::ProbabilityType::Pion,
                         achilles::Cascade::ProbabilityType::Cylinder);

    achilles::InteractionHandler interaction;
    MockEvent event;
    MockNucleus nucleus;

    SECTION("Must have exactly one kicked") {
        REQUIRE_CALL(event, Hadrons()).TIMES(AT_LEAST(2)).LR_RETURN((hadrons));

        achilles::Cascade cascade(std::move(interaction), mode, achilles::Cascade::Algorithm::Base,
                                  achilles::Cascade::InMedium::None);
        CHECK_THROWS_WITH(cascade.MeanFreePath(event, &nucleus),
                          "MeanFreePath: only one particle should be kicked.");

        cascade.SetKicked(0);
        cascade.SetKicked(1);
        CHECK_THROWS_WITH(cascade.MeanFreePath(event, &nucleus),
                          "MeanFreePath: only one particle should be kicked.");
    }

    SECTION("Must have internal test particle") {
        REQUIRE_CALL(event, Hadrons()).TIMES(AT_LEAST(1)).LR_RETURN((hadrons));
        REQUIRE_CALL(nucleus, GetPotential()).TIMES(1).RETURN(nullptr);

        achilles::Cascade cascade(std::move(interaction), mode, achilles::Cascade::Algorithm::Base,
                                  achilles::Cascade::InMedium::None);
        cascade.SetKicked(1);
        CHECK_THROWS_WITH(
            cascade.MeanFreePath(event, &nucleus),
            "MeanFreePath: kickNuc must have status -3 in order to accumulate DistanceTraveled.");
    }

    SECTION("Particle escapes marked correctly") {
        REQUIRE_CALL(event, Hadrons()).TIMES(AT_LEAST(2)).LR_RETURN((hadrons));
        REQUIRE_CALL(nucleus, Radius()).TIMES(AT_LEAST(1)).RETURN(radius);
        REQUIRE_CALL(nucleus, GetPotential()).TIMES(1).RETURN(nullptr);

        achilles::Cascade cascade(std::move(interaction), mode, achilles::Cascade::Algorithm::Base,
                                  achilles::Cascade::InMedium::None);
        cascade.SetKicked(0);
        CHECK_NOTHROW(cascade.MeanFreePath(event, &nucleus));
        CHECK(hadrons[0].Status() == achilles::ParticleStatus::final_state);
        CHECK(hadrons[0].Radius() > radius);
    }
}

/*
TEST_CASE("NuWro Mean Free Path Mode", "[Cascade]") {}

TEST_CASE("Cascade YAML", "[Cascade]") {
    auto prob = GENERATE(values<std::string>({"Gaussian", "Pion", "Cylinder"}));
    auto in_medium = GENERATE(values<std::string>({"None", "NonRelativistic", "Relativistic"}));
    YAML::Node node = YAML::Load(fmt::format(R"node(
    Interaction:
        Name: ConstantInteractions
        CrossSection: 10
    Probability: {}
    InMedium: {}
    PotentialProp: False
    Step: 0.04
    )node",
                                             prob, in_medium));

    auto cascade = node.as<achilles::Cascade>();

    CHECK(cascade.InteractionModel() == "ConstantInteractions");
    CHECK(cascade.ProbabilityModel() == prob);
    CHECK(cascade.InMediumSetting() == in_medium);
    CHECK(cascade.UsePotentialProp() == false);
    CHECK(cascade.StepSize() == 0.04);
}
*/
