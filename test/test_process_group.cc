#include "catch2/catch.hpp"

#include "Achilles/ParticleInfo.hh"
#include "Achilles/Process.hh"
#include "Achilles/XSecBackend.hh"
#include "Approx.hh"
#include "mock_classes.hh"
#include "trompeloeil.hpp"

using achilles::RegistrableBackend;

std::unique_ptr<MockBackend> MockBackend::self = nullptr;

TEST_CASE("Process Grouping Setup", "[Process]") {
    auto beam = std::make_shared<MockBeam>();
    auto nucleus = std::make_shared<MockNucleus>();

    SECTION("Setup Leptons") {
        achilles::ProcessInfo info;
        info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
        info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};

        YAML::Node config;
        auto unweight = std::make_unique<achilles::NoUnweighter>(config);
        achilles::Process process(info, std::move(unweight));
        achilles::ProcessGroup group(beam, nucleus);
        group.AddProcess(std::move(process));
        const std::vector<achilles::FourVector> momentum = {
            {1.3e+03, 0.0, 0.0, 1.3e+03},
            {1.1188e+04, 0.0, 0.0, 0.0},
            {1.27035325e+03, 6.15441682e+02, -4.52084137e+02, 1.01520877e+03},
            {1.12176467e+04, -6.15441682e+02, 4.52084137e+02, 2.84791227e+02}};

        std::vector<achilles::Particle> leptons;
        std::vector<achilles::Particle> expected_leptons;
        expected_leptons.emplace_back(info.m_leptonic.first, momentum[0]);
        expected_leptons.back().Status() = achilles::ParticleStatus::initial_state;
        expected_leptons.emplace_back(info.m_leptonic.second[0], momentum[2]);
        expected_leptons.back().Status() = achilles::ParticleStatus::final_state;

        MockEvent event;
        const MockEvent &cevent = event;
        REQUIRE_CALL(cevent, Momentum()).TIMES(1).LR_RETURN(momentum);
        REQUIRE_CALL(event, Leptons()).TIMES(1).LR_RETURN(leptons);

        group.SetupLeptons(event, std::optional<size_t>(0));

        CHECK(leptons.size() == expected_leptons.size());
        CHECK(leptons[0] == expected_leptons[0]);
        CHECK(leptons[1] == expected_leptons[1]);
    }

    SECTION("Construct Groups") {
        YAML::Node config = YAML::Load(R"config(
        Processes:
            - Leptons: [11, [11]]
        Options:
            Unweighting:
                Name: Percentile
                percentile: 99
        )config");

        MockNuclearModel model;
        achilles::ProcessInfo process_info;
        process_info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
        std::vector<achilles::ProcessInfo> infos;
        infos.push_back(process_info);
        infos.back().m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
        infos.push_back(process_info);
        infos.back().m_hadronic = {{achilles::PID::neutron()}, {achilles::PID::neutron()}};
        infos.push_back(process_info);
        infos.back().m_hadronic = {{achilles::PID::proton()},
                                   {achilles::PID::neutron(), achilles::PID::pionp()}};
        REQUIRE_CALL(model, AllowedStates(process_info)).TIMES(1).LR_RETURN(infos);
        REQUIRE_CALL(model, Mode()).TIMES(2).RETURN(achilles::NuclearMode::Quasielastic);
        REQUIRE_CALL(model, Mode()).TIMES(1).RETURN(achilles::NuclearMode::Resonance);
        auto groups = achilles::ProcessGroup::ConstructGroups(config, &model, beam, nucleus);
        CHECK(groups.size() == 2);
        CHECK(groups.at(4).Processes().size() == 2);
        CHECK(groups.at(5).Processes().size() == 1);
    }
}

TEST_CASE("Process Grouping CrossSection", "[Process]") {
    auto beam = std::make_shared<MockBeam>();
    auto nucleus = std::make_shared<MockNucleus>();
    achilles::ProcessInfo info;
    info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
    YAML::Node config;
    auto unweight = std::make_unique<achilles::NoUnweighter>(config);
    achilles::Process process(info, std::move(unweight));
    achilles::ProcessGroup group(beam, nucleus);
    group.AddProcess(std::move(process));

    YAML::Node backend_node = YAML::Load(R"backend(
    Backend:
        Name: Mock
        Options: 
    )backend");

    // MockSherpaInterface sherpa;
    auto backend = std::make_unique<MockBackend>();
    auto model = std::make_unique<MockNuclearModel>();
    trompeloeil::sequence seq;
    REQUIRE_CALL(*backend, SetOptions(backend_node["Backend"]["Options"]))
        .TIMES(1)
        .IN_SEQUENCE(seq);
    REQUIRE_CALL(*backend, AddNuclearModel(trompeloeil::ne(nullptr))).TIMES(1).IN_SEQUENCE(seq);
    REQUIRE_CALL(*backend, SetSherpa(nullptr)).TIMES(1).IN_SEQUENCE(seq);
    REQUIRE_CALL(*backend, Validate()).TIMES(1).IN_SEQUENCE(seq).RETURN(true);
    REQUIRE_CALL(*backend, AddProcess(trompeloeil::_)).TIMES(1).IN_SEQUENCE(seq);

    SECTION("Cross Section") {
        MockEvent event;
        const MockEvent &cevent = event;
        double weight;
        static constexpr double expected_weight = 10;
        REQUIRE_CALL(event, Weight()).TIMES(1).LR_RETURN(weight);
        REQUIRE_CALL(*backend, CrossSection(trompeloeil::_, trompeloeil::_))
            .LR_WITH(_1 == cevent)
            .TIMES(1)
            .RETURN(expected_weight);
        MockBackend::SetSelf(std::move(backend));
        group.SetupBackend(backend_node, std::move(model), nullptr);

        SECTION("Optimize") {
            group.CrossSection(event, std::optional<size_t>());
            CHECK(weight == expected_weight);
        }

        SECTION("Generate") {
            group.CrossSection(event, std::optional<size_t>(0));
            CHECK(weight == expected_weight);
        }
    }
}

TEST_CASE("Process Grouping Single Event", "[Process]") {
    const std::vector<achilles::FourVector> momentum = {
        {1.3e+03, 0.0, 0.0, 1.3e+03},
        {1.1188e+04, 0.0, 0.0, 0.0},
        {1.27035325e+03, 6.15441682e+02, -4.52084137e+02, 1.01520877e+03},
        {1.12176467e+04, -6.15441682e+02, 4.52084137e+02, 2.84791227e+02}};
    constexpr double ps_wgt = 1;
    constexpr double flux = 1;

    std::vector<achilles::Particle> nucleons = {{achilles::PID::proton()},
                                                {achilles::PID::neutron()}};
    auto beam = std::make_shared<MockBeam>();
    auto nucleus = std::make_shared<MockNucleus>();
    REQUIRE_CALL(*nucleus, GenerateConfig()).TIMES(1);
    REQUIRE_CALL(*nucleus, Nucleons()).TIMES(AT_LEAST(1)).LR_RETURN(nucleons);
    achilles::ProcessInfo info;
    info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
    YAML::Node backend_node = YAML::Load(R"backend(
        Backend:
            Name: Mock
            Options: 
        )backend");

    // MockSherpaInterface sherpa;
    auto backend = std::make_unique<MockBackend>();
    auto model = std::make_unique<MockNuclearModel>();
    trompeloeil::sequence seq;
    REQUIRE_CALL(*backend, SetOptions(backend_node["Backend"]["Options"]))
        .TIMES(1)
        .IN_SEQUENCE(seq);
    REQUIRE_CALL(*backend, AddNuclearModel(trompeloeil::ne(nullptr))).TIMES(1).IN_SEQUENCE(seq);
    REQUIRE_CALL(*backend, SetSherpa(nullptr)).TIMES(1).IN_SEQUENCE(seq);
    REQUIRE_CALL(*backend, Validate()).TIMES(1).IN_SEQUENCE(seq).RETURN(true);
    REQUIRE_CALL(*backend, AddProcess(trompeloeil::_)).TIMES(1).IN_SEQUENCE(seq);

    SECTION("Optimize") {
        auto optimize = true;
        YAML::Node config;
        auto unweight = std::make_unique<achilles::NoUnweighter>(config);
        achilles::Process process(info, std::move(unweight));
        achilles::ProcessGroup group(beam, nucleus);
        group.AddProcess(std::move(process));

        group.SetOptimize(optimize);

        static constexpr double expected_weight = 10;
        REQUIRE_CALL(*backend, CrossSection(trompeloeil::_, trompeloeil::_))
            .TIMES(1)
            .RETURN(expected_weight);
        MockBackend::SetSelf(std::move(backend));
        group.SetupBackend(backend_node, std::move(model), nullptr);
        auto event = group.SingleEvent(momentum, ps_wgt);
        CHECK(event.Weight() == expected_weight);
    }

    SECTION("Generate") {
        auto optimize = false;
        REQUIRE_CALL(*nucleus, ProtonsIDs()).TIMES(1).RETURN(std::vector<size_t>{0});
        REQUIRE_CALL(*nucleus, NeutronsIDs()).TIMES(1).RETURN(std::vector<size_t>{1});
        REQUIRE_CALL(*nucleus, Nucleons()).TIMES(3).LR_RETURN(std::ref(nucleons));
        REQUIRE_CALL(*beam, EvaluateFlux(achilles::PID::electron(), momentum[0]))
            .TIMES(1)
            .LR_RETURN(flux);
        YAML::Node config;
        auto unweight = std::make_unique<achilles::NoUnweighter>(config);
        achilles::Process process(info, std::move(unweight));
        achilles::ProcessGroup group(beam, nucleus);
        group.AddProcess(std::move(process));

        group.SetOptimize(optimize);

        static constexpr double expected_weight = 10;
        group.MaxWeight() = 1; // Hack to ensure weight is not rescaled by 0
        REQUIRE_CALL(*backend, CrossSection(trompeloeil::_, trompeloeil::_))
            .TIMES(1)
            .RETURN(expected_weight);
        MockBackend::SetSelf(std::move(backend));
        group.SetupBackend(backend_node, std::move(model), nullptr);
        auto event = group.SingleEvent(momentum, ps_wgt);
        CHECK(event.Weight() == expected_weight);
    }
}
