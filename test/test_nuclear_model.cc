#include "Achilles/Current.hh"
#include "catch2/catch.hpp"

#include "Achilles/NuclearModel.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Units.hh"

#include "Approx.hh"
#include "mock_classes.hh"
#include "yaml-cpp/yaml.h"

using achilles::operator""_GeV;
using achilles::ProcessInfo;

TEST_CASE("CoherentModel", "[NuclearModel]") {
    YAML::Node config = YAML::Load(R"config(
NuclearModel:
    Nucleus: 1000060120
)config");
    YAML::Node ff = YAML::Load(R"form(
vector: dummy
axial: dummy
coherent: dummy
resonancevector: dummy
resonanceaxial: dummy
mecvector: dummy
mecaxial: dummy
hyperon: dummy

dummy: dummy2
)form");

    auto form_factor = std::make_unique<MockFormFactor>();

    MockFormFactorBuilder builder;
    trompeloeil::sequence seq;
    REQUIRE_CALL(builder, Vector(ff["vector"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, AxialVector(ff["axial"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, Coherent(ff["coherent"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, ResonanceVector(ff["resonancevector"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, ResonanceAxial(ff["resonanceaxial"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, MesonExchangeVector(ff["mecvector"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, MesonExchangeAxial(ff["mecaxial"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, Hyperon(ff["hyperon"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));

    auto nucleus = std::make_shared<MockNucleus>();

    SECTION("Allowed States are valid") {
        // Require to build here or else form_factor is moved before expectations are set in next
        // test
        REQUIRE_CALL(builder, build()).TIMES(1).IN_SEQUENCE(seq).LR_RETURN(std::move(form_factor));
        REQUIRE_CALL(builder, Reset()).TIMES(1).IN_SEQUENCE(seq);
        achilles::Coherent model(config, ff, builder);

        achilles::ProcessInfo info;
        info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
        auto groups = model.AllowedStates(info);
        CHECK(groups.size() == 1);
        info = groups[0];
        CHECK(info.m_hadronic ==
              ProcessInfo::hadronic_state{{achilles::PID::carbon()}, {achilles::PID::carbon()}});

        achilles::ProcessInfo invalid;
        invalid.m_leptonic = {achilles::PID::electron(), {achilles::PID::nu_electron()}};
        CHECK_THROWS_WITH(model.AllowedStates(invalid),
                          "Coherent: Requires charge 0, but found charge 1");
    }

    SECTION("CalcCurrents") {
        std::vector<achilles::FourVector> momentum = {
            {1.3e+03, 0.0, 0.0, 1.3e+03},
            {1.1188e+04, 0.0, 0.0, 0.0},
            {1.27035325e+03, 6.15441682e+02, -4.52084137e+02, 1.01520877e+03},
            {1.12176467e+04, -6.15441682e+02, 4.52084137e+02, 2.84791227e+02}};

        double Q2 = -(momentum[0] - momentum[2]).M2();
        achilles::FormFactor::Values value;
        value.Fcoh = 1;
        REQUIRE_CALL(*form_factor, call_op(Q2)).TIMES(1).LR_RETURN((value));

        // Require to build here or else form_factor is moved before expectations are set
        REQUIRE_CALL(builder, build()).TIMES(1).IN_SEQUENCE(seq).LR_RETURN(std::move(form_factor));
        REQUIRE_CALL(builder, Reset()).TIMES(1).IN_SEQUENCE(seq);
        achilles::Coherent model(config, ff, builder);

        std::vector<achilles::NuclearModel::FFInfoMap> info_map(3);
        info_map[2][achilles::PID::photon()] = {
            achilles::FormFactorInfo{achilles::FormFactorInfo::Type::FCoh, 1}};
        auto output = model.CalcCurrents({{achilles::PID::carbon(), momentum[1]}},
                                         {{achilles::PID::carbon(), momentum[3]}}, {},
                                         momentum[0] - momentum[2], info_map[2]);
        CHECK(output.size() == 1);
        auto results = output[achilles::PID::photon()];
        achilles::VCurrent expected;
        expected[0] = {momentum[1][0] + momentum[3][0], 0.0};
        expected[1] = {momentum[1][1] + momentum[3][1], 0.0};
        expected[2] = {momentum[1][2] + momentum[3][2], 0.0};
        expected[3] = {momentum[1][3] + momentum[3][3], 0.0};
        CHECK(results.size() == model.NSpins());
        CHECK(results[0].size() == 4);
        CHECK_THAT(results[0][0].real(), Catch::Matchers::WithinAbs(expected[0].real(), 1e-8));
        CHECK_THAT(results[0][0].imag(), Catch::Matchers::WithinAbs(expected[0].imag(), 1e-8));
        CHECK_THAT(results[0][1].real(), Catch::Matchers::WithinAbs(expected[1].real(), 1e-8));
        CHECK_THAT(results[0][1].imag(), Catch::Matchers::WithinAbs(expected[1].imag(), 1e-8));
        CHECK_THAT(results[0][2].real(), Catch::Matchers::WithinAbs(expected[2].real(), 1e-8));
        CHECK_THAT(results[0][2].imag(), Catch::Matchers::WithinAbs(expected[2].imag(), 1e-8));
        CHECK_THAT(results[0][3].real(), Catch::Matchers::WithinAbs(expected[3].real(), 1e-8));
        CHECK_THAT(results[0][3].imag(), Catch::Matchers::WithinAbs(expected[3].imag(), 1e-8));
    }
}

TEST_CASE("QESpectralModel", "[NuclearModel]") {
    YAML::Node config = YAML::Load(R"config(
NuclearModel:
  SpectralP: data/pke12_tot.data
  SpectralN: data/pke12_tot.data
  Ward: Coulomb
)config");
    YAML::Node ff = YAML::Load(R"form(
vector: dummy
axial: dummy
coherent: dummy
resonancevector: dummy
resonanceaxial: dummy
mecvector: dummy
mecaxial: dummy
hyperon: dummy

dummy: dummy2
)form");

    auto form_factor = std::make_unique<MockFormFactor>();

    MockFormFactorBuilder builder;
    trompeloeil::sequence seq;
    REQUIRE_CALL(builder, Vector(ff["vector"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, AxialVector(ff["axial"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, Coherent(ff["coherent"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, ResonanceVector(ff["resonancevector"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, ResonanceAxial(ff["resonanceaxial"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, MesonExchangeVector(ff["mecvector"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, MesonExchangeAxial(ff["mecaxial"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, Hyperon(ff["hyperon"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));

    auto nucleus = std::make_shared<MockNucleus>();

    SECTION("Allowed States are valid") {
        // Require to build here or else form_factor is moved before expectations are set in next
        // test
        REQUIRE_CALL(builder, build()).TIMES(1).IN_SEQUENCE(seq).LR_RETURN(std::move(form_factor));
        REQUIRE_CALL(builder, Reset()).TIMES(1).IN_SEQUENCE(seq);
        achilles::QESpectral model(config, ff, builder);

        SECTION("Neutral Current") {
            achilles::ProcessInfo info;
            info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
            // TODO: Process group return
            auto groups = model.AllowedStates(info);
            CHECK(groups.size() == 2);
            info = groups[0];
            CHECK(info.m_hadronic == ProcessInfo::hadronic_state{{achilles::PID::neutron()},
                                                                 {achilles::PID::neutron()}});
            info = groups[1];
            CHECK(info.m_hadronic == ProcessInfo::hadronic_state{{achilles::PID::proton()},
                                                                 {achilles::PID::proton()}});
        }

        SECTION("-1 Nuclear Charge") {
            achilles::ProcessInfo info;
            info.m_leptonic = {-achilles::PID::nu_electron(), {-achilles::PID::electron()}};
            // TODO: Process group return
            // model.AllowedStates(info);
            // CHECK(info.m_states[{achilles::PID::proton()}] ==
            // std::vector<achilles::PID>{achilles::PID::neutron()});
        }

        SECTION("+1 Nuclear Charge") {
            achilles::ProcessInfo info;
            info.m_leptonic = {achilles::PID::nu_electron(), {achilles::PID::electron()}};
            // TODO: Process group return
            // model.AllowedStates(info);
            // CHECK(info.m_states[{achilles::PID::neutron()}] ==
            // std::vector<achilles::PID>{achilles::PID::proton()});
        }

        SECTION("-2 Nuclear Charge") {
            achilles::ProcessInfo invalid;
            invalid.m_leptonic = {achilles::PID::electron(), {-achilles::PID::electron()}};
            CHECK_THROWS_WITH(model.AllowedStates(invalid),
                              "Quasielastic: Requires |charge| < 2, but found |charge| 2");
        }

        SECTION("+2 Nuclear Charge") {
            achilles::ProcessInfo invalid;
            invalid.m_leptonic = {-achilles::PID::electron(), {achilles::PID::electron()}};
            CHECK_THROWS_WITH(model.AllowedStates(invalid),
                              "Quasielastic: Requires |charge| < 2, but found |charge| 2");
        }
    }

    // TODO: This test is too sensitive to minor numerical changes. Need to find a way to better
    // test
    SECTION("CalcCurrents") {
        std::vector<achilles::FourVector> momentum = {
            {1.30000000e+03, 0.00000000e+00, 0.00000000e+00, 1.30000000e+03},
            {6.61445463e+02, 5.72268451e+01, -5.50505424e+02, 1.44888604e+02},
            {6.05909852e+02, -4.66914933e+01, -3.46184656e+02, 4.95078618e+02},
            {1.35553561e+03, 1.03918338e+02, -2.04320768e+02, 9.49809986e+02}};
        double Q2 = -(momentum[0] - momentum[2]).M2() / 1.0_GeV / 1.0_GeV;
        achilles::FormFactor::Values value;
        value.F1p = 1;
        value.F1n = 1;
        value.F2p = 1;
        value.F2n = 1;
        value.FA = 1;
        REQUIRE_CALL(*form_factor, call_op(Q2)).TIMES(AT_LEAST(1)).LR_RETURN((value));

        // Require to build here or else form_factor is moved before expectations are set in next
        // test
        REQUIRE_CALL(builder, build()).TIMES(1).IN_SEQUENCE(seq).LR_RETURN(std::move(form_factor));
        REQUIRE_CALL(builder, Reset()).TIMES(1).IN_SEQUENCE(seq);
        achilles::QESpectral model(config, ff, builder);

        std::vector<achilles::NuclearModel::FFInfoMap> info_map(3);
        info_map[0][achilles::PID::photon()] = {
            achilles::FormFactorInfo{achilles::FormFactorInfo::Type::F1p, 1},
            achilles::FormFactorInfo{achilles::FormFactorInfo::Type::F2p, 1},
            achilles::FormFactorInfo{achilles::FormFactorInfo::Type::FA, 1}};
        info_map[1][achilles::PID::photon()] = {
            achilles::FormFactorInfo{achilles::FormFactorInfo::Type::F1n, 1},
            achilles::FormFactorInfo{achilles::FormFactorInfo::Type::F2n, 1},
            achilles::FormFactorInfo{achilles::FormFactorInfo::Type::FA, 1}};
        auto output_p = model.CalcCurrents({{2212, momentum[1]}}, {{2212, momentum[3]}}, {},
                                           momentum[0] - momentum[2], info_map[0]);
        auto output_n = model.CalcCurrents({{2112, momentum[1]}}, {{2112, momentum[3]}}, {},
                                           momentum[0] - momentum[2], info_map[1]);
        CHECK(output_p.size() == 1);
        auto results_p = output_p[achilles::PID::photon()];
        CHECK(results_p.size() == model.NSpins());
        CHECK(results_p[0].size() == 4);
        CHECK(output_n.size() == 1);
        auto results_n = output_n[achilles::PID::photon()];
        CHECK(results_n.size() == model.NSpins());
        CHECK(results_n[0].size() == 4);

        // TODO: Add more robust tests
    }
}
