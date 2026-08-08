#include "Achilles/Settings.hh"
#include "catch2/catch.hpp"

#include "Achilles/ParticleInfo.hh"
#include "Achilles/Process.hh"
#include "Achilles/Random.hh"
#include "Achilles/XSecBackend.hh"
#include "mock_classes.hh"

#include <array>
#include <numeric>

using achilles::Particle;
using achilles::RegistrableBackend;

std::unique_ptr<MockBackend> MockBackend::self = nullptr;

// A process whose only relevant property here is its incoming lepton species.
static achilles::Process MakeProcess(achilles::PID nu) {
    achilles::ProcessInfo info;
    info.m_leptonic = {nu, {nu}};
    YAML::Node config;
    return achilles::Process(info, std::make_unique<achilles::NoUnweighter>(config));
}

TEST_CASE("Geometry process selection", "[Process][Geometry]") {
    auto beam = std::make_shared<achilles::Beam>();
    auto nucleus = MakeNucleus();

    // Group with three processes: numu, numu, nue.
    achilles::ProcessGroup group(beam, nucleus);
    group.AddProcess(MakeProcess(achilles::PID(14)));
    group.AddProcess(MakeProcess(achilles::PID(14)));
    group.AddProcess(MakeProcess(achilles::PID(12)));
    group.ProcessWeights() = {0.5, 0.3, 0.2}; // normalized within the group
    group.MaxWeight() = 10.0;

    SECTION("NeutrinoMaxWeight sums only the matching species") {
        CHECK(group.NeutrinoMaxWeight(achilles::PID(14)) == Approx(8.0));
        CHECK(group.NeutrinoMaxWeight(achilles::PID(12)) == Approx(2.0));
        CHECK(group.NeutrinoMaxWeight(achilles::PID(16)) == 0.0);
    }

    SECTION("SelectProcess restricts to the requested species, weighted") {
        achilles::Random::Instance().Seed(12345);
        std::array<int, 3> counts{0, 0, 0};
        const int N = 400000;
        for(int i = 0; i < N; ++i) counts[group.SelectProcess(achilles::PID(14))]++;

        CHECK(counts[2] == 0); // nue process never chosen for a numu ray
        // The two numu processes are chosen in proportion to their weights.
        CHECK(static_cast<double>(counts[0]) / counts[1] == Approx(0.5 / 0.3).epsilon(0.02));

        // A nue ray always selects the single nue process.
        for(int i = 0; i < 100; ++i) CHECK(group.SelectProcess(achilles::PID(12)) == 2);
    }
}

TEST_CASE("Geometry group selection", "[Process][Geometry]") {
    auto beam = std::make_shared<achilles::Beam>();
    const achilles::PID carbon{1000060120}, hydrogen{1000010010};

    // Group A (carbon): tiny numu weight but huge nue weight.
    auto nucA = MakeNucleus(6, 12);
    achilles::ProcessGroup groupA(beam, nucA);
    groupA.AddProcess(MakeProcess(achilles::PID(14)));
    groupA.AddProcess(MakeProcess(achilles::PID(12)));
    groupA.ProcessWeights() = {0.01, 0.99};
    groupA.MaxWeight() = 100.0; // numu maxweight = 1.0, nue maxweight = 99

    // Group B (carbon): numu only.
    auto nucB = MakeNucleus(6, 12);
    achilles::ProcessGroup groupB(beam, nucB);
    groupB.AddProcess(MakeProcess(achilles::PID(14)));
    groupB.ProcessWeights() = {1.0};
    groupB.MaxWeight() = 50.0; // numu maxweight = 50

    // Group C (hydrogen): numu, but the wrong target.
    auto nucC = MakeNucleus(1, 1);
    achilles::ProcessGroup groupC(beam, nucC);
    groupC.AddProcess(MakeProcess(achilles::PID(14)));
    groupC.ProcessWeights() = {1.0};
    groupC.MaxWeight() = 1000.0;

    std::vector<achilles::ProcessGroup> groups;
    groups.push_back(std::move(groupA));
    groups.push_back(std::move(groupB));
    groups.push_back(std::move(groupC));

    SECTION("Weights use only the matching species on the target nucleus") {
        auto w = achilles::GeometryGroupWeights(groups, achilles::PID(14), carbon);
        REQUIRE(w.size() == 3);
        CHECK(w[0] == Approx(1.0));  // A: only the numu part, not its total of 100
        CHECK(w[1] == Approx(50.0)); // B
        CHECK(w[2] == 0.0);          // C: wrong target
        // The numu-on-carbon envelope (what EventGen::VertexEnvelope sums).
        CHECK(std::accumulate(w.begin(), w.end(), 0.0) == Approx(51.0));
    }

    SECTION("Selection frequencies follow the matching weights") {
        auto w = achilles::GeometryGroupWeights(groups, achilles::PID(14), carbon);
        achilles::Random::Instance().Seed(987);
        std::array<int, 3> counts{0, 0, 0};
        const int N = 400000;
        for(int i = 0; i < N; ++i) counts[achilles::Random::Instance().SelectIndex(w)]++;

        CHECK(counts[2] == 0); // wrong-target group never selected
        // B is chosen ~50x more than A (50.0 : 1.0), despite A's large nue weight.
        CHECK(static_cast<double>(counts[1]) / counts[0] == Approx(50.0).epsilon(0.05));
    }
}

TEST_CASE("Process Grouping Setup", "[Process]") {
    auto beam = std::make_shared<achilles::Beam>();
    auto nucleus = MakeNucleus();

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

        achilles::Event event;
        event.Momentum() = momentum;
        event.Leptons() = leptons;

        group.SetupLeptons(event, std::optional<size_t>(0));

        leptons = event.Leptons();
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
        auto groups = achilles::ProcessGroup::ConstructGroups(achilles::Settings(config), &model,
                                                              beam, nucleus);
        CHECK(groups.size() == 2);
        CHECK(groups.at(4).Processes().size() == 2);
        CHECK(groups.at(5).Processes().size() == 1);
    }
}

TEST_CASE("Process Grouping CrossSection", "[Process]") {
    auto beam = std::make_shared<achilles::Beam>();
    auto nucleus = MakeNucleus();
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
        achilles::Event event;
        static constexpr double expected_weight = 10;
        REQUIRE_CALL(*backend, CrossSection(trompeloeil::_, trompeloeil::_))
            .LR_WITH(_1 == event)
            .TIMES(1)
            .RETURN(expected_weight);
        MockBackend::SetSelf(std::move(backend));
        group.SetupBackend(backend_node, std::move(model), nullptr);

        SECTION("Optimize") {
            group.CrossSection(event, std::optional<size_t>());
            CHECK(event.Weight() == expected_weight);
        }

        SECTION("Generate") {
            group.CrossSection(event, std::optional<size_t>(0));
            CHECK(event.Weight() == expected_weight);
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

    std::vector<achilles::Particle> nucleons = {Particle{achilles::PID::proton()},
                                                Particle{achilles::PID::neutron()}};
    auto mock_flux = std::make_shared<MockFluxType>();
    auto beam = MakeBeam(achilles::PID::electron(), mock_flux);
    auto density = std::make_unique<MockDensity>();
    REQUIRE_CALL(*density, GetConfiguration()).TIMES(1).LR_RETURN(nucleons);
    auto nucleus = MakeNucleus();
    nucleus->SetDensity(std::move(density));

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
        REQUIRE_CALL(*mock_flux, EvaluateFlux(momentum[0]))
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
