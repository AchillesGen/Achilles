#include "catch2/catch.hpp"

#include "mock_classes.hh"
#include <limits>

#include "Achilles/Cascade.hh"
#include "Achilles/CascadeInteractions/ConstantInteractions.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventHistory.hh"
#include "Achilles/FoldingCorrection.hh"
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

namespace {

// Elastic NN scattering with a constant 10 mb cross section, i.e. sigma = 1 fm^2
achilles::InteractionHandler ConstantNN() {
    achilles::InteractionHandler handler;
    achilles::ConstantInteraction interaction;
    const auto p = achilles::PID::proton(), n = achilles::PID::neutron();
    interaction.AddInteraction({p, p}, {{{p, p}, 10}});
    interaction.AddInteraction({p, n}, {{{p, n}, 10}});
    interaction.AddInteraction({n, n}, {{{n, n}, 10}});
    handler.RegisterInteraction(std::make_unique<achilles::ConstantInteraction>(interaction));
    return handler;
}

} // namespace

TEST_CASE("Veto cascade", "[Cascade]") {
    const achilles::Particle projectile{achilles::PID::proton(),
                                        {1000, 500, 0, 0},
                                        {0, 0, 0},
                                        achilles::ParticleStatus::propagating};
    constexpr double radius = 5;

    SECTION("Escapes without targets") {
        achilles::Particles hadrons = {projectile};
        MockEvent event;
        MockNucleus nucleus;
        ALLOW_CALL(event, Hadrons()).LR_RETURN((hadrons));
        ALLOW_CALL(nucleus, Radius()).RETURN(radius);

        achilles::Cascade cascade(ConstantNN(), achilles::Cascade::ProbabilityType::Gaussian,
                                  achilles::Cascade::Algorithm::Veto,
                                  achilles::Cascade::InMedium::None);
        cascade.Evolve(event, &nucleus);

        CHECK(hadrons[0].Status() == achilles::ParticleStatus::final_state);
        CHECK(hadrons[0].Radius() > radius);
    }

    SECTION("No collision inside the formation zone") {
        achilles::Particles hadrons = {projectile, Particle{achilles::PID::neutron(),
                                                            {achilles::Constant::mN, 0, 0, 0},
                                                            {1, 0, 0},
                                                            achilles::ParticleStatus::background}};
        hadrons[0].SetFormationZone({10000, 0, 0, 0}, {88.2, 0, 0, 0});
        MockEvent event;
        MockNucleus nucleus;
        ALLOW_CALL(event, Hadrons()).LR_RETURN((hadrons));
        ALLOW_CALL(nucleus, Radius()).RETURN(radius);

        achilles::Cascade cascade(ConstantNN(), achilles::Cascade::ProbabilityType::Cylinder,
                                  achilles::Cascade::Algorithm::Veto,
                                  achilles::Cascade::InMedium::None);
        cascade.Evolve(event, &nucleus);

        CHECK(hadrons[0].Status() == achilles::ParticleStatus::final_state);
        CHECK(hadrons[1].Status() == achilles::ParticleStatus::background);
    }

    SECTION("Potential propagation is not supported") {
        achilles::Particles hadrons = {projectile};
        MockEvent event;
        MockNucleus nucleus;
        ALLOW_CALL(event, Hadrons()).LR_RETURN((hadrons));
        ALLOW_CALL(nucleus, Radius()).RETURN(radius);

        achilles::Cascade cascade(ConstantNN(), achilles::Cascade::ProbabilityType::Gaussian,
                                  achilles::Cascade::Algorithm::Veto,
                                  achilles::Cascade::InMedium::None, true);
        CHECK_THROWS_AS(cascade.Evolve(event, &nucleus), std::runtime_error);
    }

    SECTION("NaN cross section never interacts") {
        // A NaN probability must count as no interaction, as in the Base algorithm
        achilles::InteractionHandler handler;
        achilles::ConstantInteraction interaction;
        interaction.AddInteraction({achilles::PID::proton(), achilles::PID::neutron()},
                                   {{{achilles::PID::proton(), achilles::PID::neutron()},
                                     std::numeric_limits<double>::quiet_NaN()}});
        handler.RegisterInteraction(std::make_unique<achilles::ConstantInteraction>(interaction));
        achilles::Particles hadrons = {projectile, Particle{achilles::PID::neutron(),
                                                            {achilles::Constant::mN, 0, 0, 0},
                                                            {2, 0.3, 0},
                                                            achilles::ParticleStatus::background}};
        MockEvent event;
        MockNucleus nucleus;
        ALLOW_CALL(event, Hadrons()).LR_RETURN((hadrons));
        ALLOW_CALL(nucleus, Radius()).RETURN(radius);
        ALLOW_CALL(nucleus, ID()).RETURN(achilles::PID::hydrogen());

        achilles::Cascade cascade(std::move(handler), achilles::Cascade::ProbabilityType::Gaussian,
                                  achilles::Cascade::Algorithm::Veto,
                                  achilles::Cascade::InMedium::None);
        CHECK_NOTHROW(cascade.Evolve(event, &nucleus));
        CHECK(hadrons[1].Status() == achilles::ParticleStatus::background);
    }

    SECTION("SurvivalProfile matches the first-collision times") {
        const std::vector<achilles::Particle> targets = {
            Particle{achilles::PID::neutron(),
                     {achilles::Constant::mN, 0, 0, 0},
                     {1, 0.2, 0},
                     achilles::ParticleStatus::background},
            Particle{achilles::PID::neutron(),
                     {achilles::Constant::mN, 0, 0, 0},
                     {2, -0.4, 0.1},
                     achilles::ParticleStatus::background},
            Particle{achilles::PID::neutron(),
                     {achilles::Constant::mN, 0, 0, 0},
                     {3, 0.1, -0.3},
                     achilles::ParticleStatus::background}};
        achilles::Particles hadrons;
        MockEvent event;
        MockNucleus nucleus;
        ALLOW_CALL(event, Hadrons()).LR_RETURN((hadrons));
        ALLOW_CALL(nucleus, Radius()).RETURN(radius);
        ALLOW_CALL(nucleus, ID()).RETURN(achilles::PID::hydrogen());

        achilles::Cascade cascade(ConstantNN(), achilles::Cascade::ProbabilityType::Gaussian,
                                  achilles::Cascade::Algorithm::Veto,
                                  achilles::Cascade::InMedium::None);
        // Projectile velocity 0.5: the targets are passed at t = 2, 4, 6 fm
        const std::vector<double> times = {1, 3, 5, 7};
        hadrons = {projectile};
        hadrons.insert(hadrons.end(), targets.begin(), targets.end());
        const auto expected = cascade.SurvivalProfile(event, &nucleus, 0, times);

        constexpr size_t ntrials = 20000;
        std::vector<double> survived(times.size(), 0);
        for(size_t i = 0; i < ntrials; ++i) {
            hadrons = {projectile};
            hadrons.insert(hadrons.end(), targets.begin(), targets.end());
            cascade.Evolve(event, &nucleus);
            const auto &first = cascade.FirstCollisions();
            const auto it = first.find(0);
            for(size_t k = 0; k < times.size(); ++k)
                survived[k] += it == first.end() || it->second.time > times[k];
        }
        CHECK(expected[0] == Approx(1));
        CHECK(expected.back() < 0.9);
        for(size_t k = 0; k < times.size(); ++k) {
            const double frac = survived[k] / ntrials;
            const double err = std::sqrt(expected[k] * (1 - expected[k]) / ntrials);
            CHECK(frac == Approx(expected[k]).margin(5 * err + 1e-9));
        }
    }

    SECTION("Single pair interacts with probability P(b^2, sigma)") {
        // Hydrogen skips Pauli blocking, so the pair's interaction probability is exactly P
        const double b = GENERATE(0.3, 0.7);
        auto mode = GENERATE(achilles::Cascade::ProbabilityType::Gaussian,
                             achilles::Cascade::ProbabilityType::Cylinder);
        auto alg = GENERATE(achilles::Cascade::Algorithm::Veto, achilles::Cascade::Algorithm::Base);
        const double sigma = 1; // fm^2
        const double expected = mode == achilles::Cascade::ProbabilityType::Gaussian
                                    ? std::exp(-M_PI * b * b / sigma)
                                    : (b * b < sigma / M_PI ? 1.0 : 0.0);

        const achilles::Particle target{achilles::PID::neutron(),
                                        {achilles::Constant::mN, 0, 0, 0},
                                        {2, b, 0},
                                        achilles::ParticleStatus::background};
        achilles::Particles hadrons;
        MockEvent event;
        MockNucleus nucleus;
        ALLOW_CALL(event, Hadrons()).LR_RETURN((hadrons));
        ALLOW_CALL(nucleus, Radius()).RETURN(radius);
        ALLOW_CALL(nucleus, ID()).RETURN(achilles::PID::hydrogen());
        ALLOW_CALL(nucleus, GetPotential()).RETURN(nullptr);

        achilles::Cascade cascade(ConstantNN(), mode, alg, achilles::Cascade::InMedium::None);
        constexpr size_t ntrials = 4000;
        size_t hits = 0;
        for(size_t i = 0; i < ntrials; ++i) {
            hadrons = {projectile, target};
            cascade.Evolve(event, &nucleus);
            hits += hadrons[1].Status() == achilles::ParticleStatus::interacted;
        }
        const double frac = static_cast<double>(hits) / ntrials;
        const double err = std::sqrt(expected * (1 - expected) / ntrials);
        CHECK(frac == Approx(expected).margin(5 * err + 1e-9));
    }
}

TEST_CASE("FoldingCorrection shifts the outgoing lepton", "[Cascade]") {
    const Particle beam{
        achilles::PID::electron(), {1000, 0, 0, 1000}, {}, achilles::ParticleStatus::initial_state};
    const Particle scattered{
        achilles::PID::electron(), {800, 0, 480, 640}, {}, achilles::ParticleStatus::final_state};
    const Particle nucleon{
        achilles::PID::proton(), {1100, 0, -480, 360}, {}, achilles::ParticleStatus::final_state};
    achilles::Event event;
    event.Leptons() = {beam, scattered};
    event.Hadrons() = {nucleon};
    event.History().AddVertex({}, {beam}, {scattered, nucleon},
                              achilles::EventHistoryNode::StatusCode::primary);

    SECTION("Energy lowered at fixed direction, in the event and the history") {
        REQUIRE(achilles::FoldingCorrection::ShiftLepton(event, 30));
        const auto &lep = event.Leptons()[1].Momentum();
        const double mass = event.Leptons()[1].Info().Mass();
        CHECK(lep.E() == Approx(770));
        CHECK(lep.P() == Approx(std::sqrt(770 * 770 - mass * mass)));
        CHECK(lep.Py() / lep.Pz() == Approx(480.0 / 640.0));
        CHECK(event.Leptons()[0].Momentum().E() == Approx(1000));
        const auto &out = event.History().Primary()->ParticlesOut();
        CHECK(out[0].Momentum().E() == Approx(770));
        CHECK(out[1].Momentum().E() == Approx(1100));
    }

    SECTION("Modes parse by name") {
        using Mode = achilles::FoldingCorrection::Mode;
        CHECK(achilles::FoldingCorrection::ParseMode("Full") == Mode::Full);
        CHECK(achilles::FoldingCorrection::ParseMode("Limited") == Mode::Limited);
        CHECK(achilles::FoldingCorrection::ParseMode("Collided") == Mode::Collided);
        CHECK_THROWS(achilles::FoldingCorrection::ParseMode("Symmetric"));
    }

    SECTION("A shift below the lepton mass is refused") {
        CHECK_FALSE(achilles::FoldingCorrection::ShiftLepton(event, 900));
        CHECK(event.Leptons()[1].Momentum().E() == Approx(800));
        CHECK(event.History().Primary()->ParticlesOut()[0].Momentum().E() == Approx(800));
    }
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
