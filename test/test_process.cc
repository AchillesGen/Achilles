#include "catch2/catch.hpp"

#include "Achilles/ParticleInfo.hh"
#include "Achilles/Process.hh"
#include "Achilles/Unweighter.hh"
#include "Approx.hh"
#include "mock_classes.hh"
#include "trompeloeil.hpp"
#include <utility>

using achilles::Particle;

template <typename Type, std::size_t Size, std::size_t... Index>
constexpr auto Iota_helper(std::index_sequence<Index...>) {
    return std::array<Type, Size>{Index...};
}

template <typename Type, std::size_t Size> constexpr auto Iota() {
    return Iota_helper<Type, Size>(std::make_index_sequence<Size>());
}

TEST_CASE("Handles weights correctly", "[Process]") {
    achilles::ProcessInfo info;
    info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};

    YAML::Node config;
    auto unweight = std::make_unique<achilles::NoUnweighter>(config);
    achilles::Process process(info, std::move(unweight));
    static constexpr size_t ncalls = 10;
    static constexpr std::array<double, ncalls> weights = Iota<double, ncalls>();

    for(const auto &weight : weights) { process.AddWeight(weight); }

    static constexpr double expected_xsec = static_cast<double>(ncalls - 1) / 2;
    CHECK_THAT(process.TotalCrossSection(), Catch::Matchers::WithinRel(expected_xsec));

    static constexpr double max_wgt = weights.back();
    CHECK_THAT(process.MaxWeight(), Catch::Matchers::WithinRel(max_wgt));

    static constexpr double value = 42;
    CHECK_THAT(process.Unweight(value), Catch::Matchers::WithinRel(value));
    CHECK_THAT(process.UnweightEff(), Catch::Matchers::WithinRel(expected_xsec / max_wgt));
}

TEST_CASE("Handles events correctly", "[Process]") {
    achilles::ProcessInfo info;
    info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
    const std::vector<achilles::FourVector> momentum = {
        {1.3e+03, 0.0, 0.0, 1.3e+03},
        {1.1188e+04, 0.0, 0.0, 0.0},
        {1.27035325e+03, 6.15441682e+02, -4.52084137e+02, 1.01520877e+03},
        {1.12176467e+04, -6.15441682e+02, 4.52084137e+02, 2.84791227e+02}};

    SECTION("Extract Momentum") {
        info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
        achilles::FourVector lep_in, lep_in_expected{1.3e+03, 0.0, 0.0, 1.3e+03};
        std::vector<achilles::FourVector> had_in, lep_out, had_out, spect;
        std::vector<achilles::FourVector> had_in_exp{{1.1188e+04, 0.0, 0.0, 0.0}};
        std::vector<achilles::FourVector> lep_out_exp{
            {1.27035325e+03, 6.15441682e+02, -4.52084137e+02, 1.01520877e+03}};
        std::vector<achilles::FourVector> had_out_exp{
            {1.12176467e+04, -6.15441682e+02, 4.52084137e+02, 2.84791227e+02}};
        YAML::Node config;
        auto unweight = std::make_unique<achilles::NoUnweighter>(config);
        achilles::Process process(info, std::move(unweight));

        const MockEvent event{};
        REQUIRE_CALL(event, Momentum()).TIMES(1).RETURN(momentum);

        process.ExtractMomentum(event, lep_in, had_in, lep_out, had_out, spect);
        CHECK(lep_in == lep_in_expected);
        CHECK_THAT(had_in, AllFourVectorApprox(had_in_exp));
        CHECK_THAT(lep_out, AllFourVectorApprox(lep_out_exp));
        CHECK_THAT(had_out, AllFourVectorApprox(had_out_exp));
    }

    SECTION("Extract Particles") {
        info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
        auto expected = momentum[0] - momentum[2];
        YAML::Node config;
        auto unweight = std::make_unique<achilles::NoUnweighter>(config);
        achilles::Process process(info, std::move(unweight));

        const MockEvent event{};
        REQUIRE_CALL(event, Momentum()).TIMES(1).RETURN(momentum);

        auto qvec = process.ExtractQ(event);
        CHECK_THAT(qvec, FourVectorApprox(expected));
    }

    SECTION("Extract Q Vector") {
        info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
        achilles::Particle lep_in;
        achilles::Particle lep_in_expected{achilles::PID::electron(), {1.3e+03, 0.0, 0.0, 1.3e+03}};
        std::vector<achilles::Particle> had_in, lep_out, had_out, spect;
        std::vector<achilles::Particle> had_in_exp{
            Particle{achilles::PID::proton(), {1.1188e+04, 0.0, 0.0, 0.0}}};
        std::vector<achilles::Particle> lep_out_exp{
            Particle{achilles::PID::electron(),
             {1.27035325e+03, 6.15441682e+02, -4.52084137e+02, 1.01520877e+03}}};
        std::vector<achilles::Particle> had_out_exp{
            Particle{achilles::PID::proton(),
             {1.12176467e+04, -6.15441682e+02, 4.52084137e+02, 2.84791227e+02}}};
        YAML::Node config;
        auto unweight = std::make_unique<achilles::NoUnweighter>(config);
        achilles::Process process(info, std::move(unweight));

        const MockEvent event{};
        REQUIRE_CALL(event, Momentum()).TIMES(1).RETURN(momentum);

        process.ExtractParticles(event, lep_in, had_in, lep_out, had_out, spect);
        CHECK(lep_in == lep_in_expected);
        CHECK_THAT(had_in, AllParticleApprox(had_in_exp));
        CHECK_THAT(lep_out, AllParticleApprox(lep_out_exp));
        CHECK_THAT(had_out, AllParticleApprox(had_out_exp));
    }

    SECTION("Setup Hadrons [Coherent]") {
        info.m_hadronic = {{achilles::PID::carbon()}, {achilles::PID::carbon()}};
        YAML::Node config;
        auto unweight = std::make_unique<achilles::NoUnweighter>(config);
        achilles::Process process(info, std::move(unweight));

        std::vector<achilles::Particle> nucleons;
        std::vector<achilles::Particle> expected_nucleons;
        expected_nucleons.emplace_back(achilles::PID::carbon(), momentum[1]);
        expected_nucleons.back().Status() = achilles::ParticleStatus::initial_state;
        expected_nucleons.emplace_back(achilles::PID::carbon(), momentum[3]);
        expected_nucleons.back().Status() = achilles::ParticleStatus::final_state;

        MockEvent event;
        const MockEvent &cevent = event;
        REQUIRE_CALL(cevent, Momentum()).TIMES(1).LR_RETURN(momentum);
        REQUIRE_CALL(event, Hadrons()).TIMES(4).LR_RETURN(std::ref(nucleons));

        process.SetupHadrons(event);

        CHECK(nucleons.size() == expected_nucleons.size());
        CHECK(nucleons[0] == expected_nucleons[0]);
        CHECK(nucleons[1] == expected_nucleons[1]);
    }

    SECTION("Setup Hadrons [Quasielastic]") {
        info.m_hadronic =
            GENERATE(as<achilles::ProcessInfo::hadronic_state>{},
                     achilles::ProcessInfo::hadronic_state{{achilles::PID::proton()},
                                                           {achilles::PID::proton()}},
                     achilles::ProcessInfo::hadronic_state{{achilles::PID::neutron()},
                                                           {achilles::PID::neutron()}},
                     achilles::ProcessInfo::hadronic_state{{achilles::PID::proton()},
                                                           {achilles::PID::neutron()}},
                     achilles::ProcessInfo::hadronic_state{{achilles::PID::neutron()},
                                                           {achilles::PID::proton()}});
        YAML::Node config;
        auto unweight = std::make_unique<achilles::NoUnweighter>(config);
        achilles::Process process(info, std::move(unweight));

        std::vector<achilles::Particle> nucleons{Particle(achilles::PID::proton()),
                                                 Particle(achilles::PID::neutron())};
        achilles::Event event;
        event.Momentum() = momentum;
        event.Hadrons() = nucleons;

        process.SetupHadrons(event);

        std::vector<achilles::Particle> expected_nucleons;
        auto pid1 = info.m_hadronic.first[0];
        auto pid2 = info.m_hadronic.second[0];
        if(pid1 == achilles::PID::proton()) {
            expected_nucleons.emplace_back(achilles::PID::proton(), momentum[1]);
            expected_nucleons.back().Status() = achilles::ParticleStatus::initial_state;
            expected_nucleons.emplace_back(achilles::PID::neutron());
        } else {
            expected_nucleons.emplace_back(achilles::PID::proton());
            expected_nucleons.emplace_back(achilles::PID::neutron(), momentum[1]);
            expected_nucleons.back().Status() = achilles::ParticleStatus::initial_state;
        }
        expected_nucleons.emplace_back(pid2, momentum[3]);
        expected_nucleons.back().Status() = achilles::ParticleStatus::propagating;

        nucleons = event.Hadrons();
        CHECK(nucleons.size() == expected_nucleons.size());
        CHECK(nucleons[0] == expected_nucleons[0]);
        CHECK(nucleons[1] == expected_nucleons[1]);
        CHECK(nucleons[2] == expected_nucleons[2]);
    }
}

TEST_CASE("Correct metadata", "[Process]") {
    achilles::ProcessInfo info;
    info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
    YAML::Node config;
    auto unweight = std::make_unique<achilles::NoUnweighter>(config);
    achilles::Process process(info, std::move(unweight));

    MockBackend backend;
    MockNuclearModel model;

    REQUIRE_CALL(model, Mode()).TIMES(1).RETURN(achilles::NuclearMode::Quasielastic);
    REQUIRE_CALL(model, GetName()).TIMES(2).RETURN("QESpectral");
    REQUIRE_CALL(model, InspireHEP()).TIMES(1).RETURN("Isaacson:2022cwh");
    REQUIRE_CALL(backend, GetNuclearModel()).TIMES(3).LR_RETURN(&model);

    process.SetID(&model);
    auto metadata = process.Metadata(&backend);

    CHECK(metadata.id == 251);
    CHECK(metadata.name == "QESpectralNC1p0pi");
    auto msg = fmt::format("QESpectral ProcessInfo([11, 2212] -> [11, 2212])");
    CHECK(metadata.description == msg);
    CHECK(metadata.inspireHEP == "Isaacson:2022cwh");
}
