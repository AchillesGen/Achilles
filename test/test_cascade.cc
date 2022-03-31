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

        nuchic::Cascade cascade(std::move(interaction), nuchic::Cascade::ProbabilityType::Gaussian,
                                nuchic::Cascade::InMedium::None);
        cascade.Kick(nucleus, {0, 100, 0, 0}, {10, 0});
        CHECK(particles[0].Status() == nuchic::ParticleStatus::propagating);
        CHECK(particles[1].Status() == nuchic::ParticleStatus::background);
        particles[0].Status() = nuchic::ParticleStatus::background;

        cascade.Reset();
        cascade.Kick(nucleus, {0, 100, 0, 0}, {0, 10});
        CHECK(particles[0].Status() == nuchic::ParticleStatus::background);
        CHECK(particles[1].Status() == nuchic::ParticleStatus::propagating);
    }
}

TEST_CASE("Evolve States: 1 nucleon", "[Cascade]") {
    nuchic::Particles hadrons = {{nuchic::PID::proton(), {1000, 100, 0, 0},
                                 {0, 0, 0}, nuchic::ParticleStatus::propagating}};
    constexpr double radius = 1;

    auto mode = GENERATE(nuchic::Cascade::ProbabilityType::Gaussian,
                         nuchic::Cascade::ProbabilityType::Pion,
                         nuchic::Cascade::ProbabilityType::Cylinder);

    SECTION("Captured") {
        auto interaction = std::make_unique<MockInteraction>();
        auto nucleus = std::make_shared<MockNucleus>();
        auto potential = std::make_shared<MockPotential>();

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));
        REQUIRE_CALL(*nucleus, GetPotential())
            .TIMES(1)
            .LR_RETURN((potential));

        REQUIRE_CALL(*potential, Hamiltonian(hadrons[0].Momentum().P(), hadrons[0].Position().P()))
            .TIMES(1)
            .RETURN(5);

        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None, true);
        cascade.SetKicked(0);
        cascade.Evolve(nucleus);

        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::captured);
    }

    SECTION("Large Formation Zone") {
        auto interaction = std::make_unique<MockInteraction>();
        auto nucleus = std::make_shared<MockNucleus>();
        auto potential = std::make_shared<MockPotential>();

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));
        REQUIRE_CALL(*nucleus, GetPotential())
            .TIMES(1)
            .RETURN(nullptr);

        REQUIRE_CALL(*nucleus, Radius())
            .TIMES(AT_LEAST(1))
            .RETURN(radius);

        hadrons[0].SetFormationZone({10000, 0, 0, 0}, {88.2, 0, 0, 0});
        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
        cascade.SetKicked(0);
        cascade.Evolve(nucleus);

        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
    }

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
        REQUIRE_CALL(*nucleus, GetPotential())
            .TIMES(1)
            .RETURN(nullptr);

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
        REQUIRE_CALL(*nucleus, GetPotential())
            .TIMES(1)
            .RETURN(nullptr);

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
        REQUIRE_CALL(*nucleus, GetPotential())
            .TIMES(1)
            .RETURN(nullptr);

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

TEST_CASE("Evolve States: 3 nucleons", "[Cascade]") {
    nuchic::Particles hadrons = {{nuchic::PID::proton(), {1000, 100, 0, 0},
                                 {0, 0, 0}, nuchic::ParticleStatus::propagating},
                                 {nuchic::PID::proton(), {nuchic::Constant::mN, 0, 0, 0},
                                 {0.5, 0, 0}, nuchic::ParticleStatus::background},
                                 {nuchic::PID::neutron(), {nuchic::Constant::mN, 0, 0, 0},
                                 {3, 0, 0}, nuchic::ParticleStatus::background}};
    constexpr double radius = 4;

    auto mode = GENERATE(nuchic::Cascade::ProbabilityType::Gaussian,
                         nuchic::Cascade::ProbabilityType::Pion,
                         nuchic::Cascade::ProbabilityType::Cylinder);

    SECTION("Captured") {
        auto interaction = std::make_unique<MockInteraction>();
        auto nucleus = std::make_shared<MockNucleus>();
        auto potential = std::make_shared<MockPotential>();

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));
        REQUIRE_CALL(*nucleus, GetPotential())
            .TIMES(1)
            .LR_RETURN((potential));

        REQUIRE_CALL(*potential, Hamiltonian(hadrons[0].Momentum().P(), hadrons[0].Position().P()))
            .TIMES(1)
            .RETURN(5);

        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None, true);
        cascade.SetKicked(0);
        cascade.Evolve(nucleus);

        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::captured);
        CHECK(hadrons[1].Status() == nuchic::ParticleStatus::background);
        CHECK(hadrons[2].Status() == nuchic::ParticleStatus::background);
    }

    SECTION("Large Formation Zone") {
        auto interaction = std::make_unique<MockInteraction>();
        auto nucleus = std::make_shared<MockNucleus>();
        auto potential = std::make_shared<MockPotential>();

        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(2)
            .LR_RETURN((hadrons));
        REQUIRE_CALL(*nucleus, GetPotential())
            .TIMES(1)
            .RETURN(nullptr);

        REQUIRE_CALL(*nucleus, Radius())
            .TIMES(AT_LEAST(1))
            .RETURN(radius);

        hadrons[0].SetFormationZone({10000, 0, 0, 0}, {88.2, 0, 0, 0});
        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
        cascade.SetKicked(0);
        cascade.Evolve(nucleus);

        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
        CHECK(hadrons[1].Status() == nuchic::ParticleStatus::background);
        CHECK(hadrons[2].Status() == nuchic::ParticleStatus::background);
    }

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
        REQUIRE_CALL(*nucleus, GetPotential())
            .TIMES(AT_LEAST(1))
            .RETURN(nullptr);
        REQUIRE_CALL(*nucleus, Rho(trompeloeil::gt(0)))
            .TIMES(4)
            .RETURN(0);

        REQUIRE_CALL(*nucleus, Radius())
            .TIMES(AT_LEAST(1))
            .RETURN(radius);

        std::pair<nuchic::FourVector, nuchic::FourVector> output{{1000, 80, 0, 0}, {150, -20, 0, 0}};
        REQUIRE_CALL(*interaction, CrossSection(trompeloeil::_, hadrons[1]))
            .TIMES(1)
            .RETURN(1000);
        REQUIRE_CALL(*interaction, CrossSection(trompeloeil::_, hadrons[2]))
            .TIMES(1)
            .RETURN(1000);
        REQUIRE_CALL(*interaction, FinalizeMomentum(trompeloeil::_, trompeloeil::_, trompeloeil::_))
            .TIMES(2)
            .LR_RETURN((output));

        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
        cascade.Evolve(&event);

        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
        CHECK(hadrons[1].Status() == nuchic::ParticleStatus::background);
        CHECK(hadrons[2].Status() == nuchic::ParticleStatus::background);
        CHECK(hadrons[0].Radius() > radius);
    }

    // SECTION("NuWro Evolve") {
    //     auto interaction = std::make_unique<MockInteraction>();
    //     auto nucleus = std::make_shared<MockNucleus>();

    //     REQUIRE_CALL(*nucleus, Nucleons())
    //         .TIMES(2)
    //         .LR_RETURN((hadrons));
    //     REQUIRE_CALL(*nucleus, GetPotential())
    //         .TIMES(1)
    //         .RETURN(nullptr);

    //     REQUIRE_CALL(*nucleus, Radius())
    //         .TIMES(AT_LEAST(1))
    //         .RETURN(radius);

    //     nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
    //     cascade.SetKicked(0);
    //     cascade.NuWro(nucleus);

    //     CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
    //     CHECK(hadrons[0].Radius() > radius);
    // }
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
        REQUIRE_CALL(*nucleus, GetPotential())
            .TIMES(1)
            .RETURN(nullptr);

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
        REQUIRE_CALL(*nucleus, GetPotential())
            .TIMES(1)
            .RETURN(nullptr);

        nuchic::Cascade cascade(std::move(interaction), mode, nuchic::Cascade::InMedium::None);
        cascade.SetKicked(0);
        CHECK_NOTHROW(cascade.MeanFreePath(nucleus));
        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::escaped);
        CHECK(hadrons[0].Radius() > radius);
    }
}

TEST_CASE("NuWro Mean Free Path Mode", "[Cascade]") {

}

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
    )node", prob, in_medium));

    auto cascade = node.as<nuchic::Cascade>();

    CHECK(cascade.InteractionModel() == "ConstantInteractions");
    CHECK(cascade.ProbabilityModel() == prob);
    CHECK(cascade.InMediumSetting() == in_medium);
    CHECK(cascade.UsePotentialProp() == false);
    CHECK(cascade.StepSize() == 0.04);
}
