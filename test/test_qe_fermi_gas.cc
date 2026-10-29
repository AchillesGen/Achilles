// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/Particle.hh"
#include "Achilles/QEFermiGas.hh"

#include "mock_classes.hh"
#include "yaml-cpp/yaml.h"

#include <cmath>
#include <random>

using Catch::Matchers::WithinRel;
using trompeloeil::_;

namespace {

const std::string kCFG = "{Distribution: CFG, kF: 220, Binding: 25, c0: 4.16, lambda: 2.75}";

std::unique_ptr<achilles::QEFermiGas> MakeModel(const std::string &fermi_gas) {
    YAML::Node config = YAML::Load("NuclearModel:\n  Ward: None\n  FermiGas: " + fermi_gas);
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
    MockFormFactorBuilder builder;
    ALLOW_CALL(builder, Vector(_, _)).LR_RETURN(std::ref(builder));
    ALLOW_CALL(builder, AxialVector(_, _)).LR_RETURN(std::ref(builder));
    ALLOW_CALL(builder, Coherent(_, _)).LR_RETURN(std::ref(builder));
    ALLOW_CALL(builder, ResonanceVector(_, _)).LR_RETURN(std::ref(builder));
    ALLOW_CALL(builder, ResonanceAxial(_, _)).LR_RETURN(std::ref(builder));
    ALLOW_CALL(builder, MesonExchangeVector(_, _)).LR_RETURN(std::ref(builder));
    ALLOW_CALL(builder, MesonExchangeAxial(_, _)).LR_RETURN(std::ref(builder));
    ALLOW_CALL(builder, Hyperon(_, _)).LR_RETURN(std::ref(builder));
    auto form_factor = std::make_unique<MockFormFactor>();
    REQUIRE_CALL(builder, build()).TIMES(1).LR_RETURN(std::move(form_factor));
    REQUIRE_CALL(builder, Reset()).TIMES(1);
    return std::make_unique<achilles::QEFermiGas>(config, ff, builder);
}

// 4 pi int_0^pmax p^2 Density(p) dp with a fine midpoint rule
double Normalization(const achilles::QEFermiGas &model, double pmax) {
    constexpr size_t n = 200000;
    const double h = pmax / n;
    double sum = 0;
    for(size_t i = 0; i < n; ++i) {
        const double p = (static_cast<double>(i) + 0.5) * h;
        sum += 4 * M_PI * p * p * model.Density(p) * h;
    }
    return sum;
}

} // namespace

TEST_CASE("QEFermiGas momentum distributions", "[QEFermiGas]") {
    SECTION("RFG is a normalized step") {
        auto model = MakeModel("{Distribution: RFG, kF: 221, Binding: 25}");
        CHECK_THAT(Normalization(*model, 1000), WithinRel(1.0, 1e-4));
        CHECK(model->Occupation(220) == 1);
        CHECK(model->Occupation(222) == 0);
        CHECK(model->PCut() == 221);
    }

    SECTION("CFG is normalized and has the expected tail") {
        const double c0 = 4.16, lambda = 2.75;
        auto model = MakeModel(kCFG);
        CHECK_THAT(Normalization(*model, 1000), WithinRel(1.0, 1e-4));
        const double tail = c0 / (M_PI * M_PI) * (1 - 1 / lambda);
        CHECK_THAT(model->Occupation(100), WithinRel(1 - tail, 1e-12));
        CHECK_THAT(Normalization(*model, 220), WithinRel(1 - tail, 1e-4));
        CHECK(model->Occupation(lambda * 220 + 1) == 0);
        CHECK(model->Occupation(221) < 1);
        CHECK_THAT(model->PCut(), WithinRel(lambda * 220, 1e-12));
    }
}

TEST_CASE("QEFermiGas initial state weight", "[QEFermiGas]") {
    auto model = MakeModel("{Distribution: RFG, kF: 221, Binding: 25}");
    const double norm = 3 / (4 * M_PI * std::pow(221.0, 3));
    achilles::Particle in(achilles::PID::proton(), {0, 0, 0, 100});
    achilles::Particle blocked(achilles::PID::proton(), {0, 0, 0, 150});
    achilles::Particle free(achilles::PID::proton(), {0, 0, 0, 400});
    CHECK_THAT(model->InitialStateWeight({in}, {}, 6, 6), WithinRel(6 * norm, 1e-12));
    CHECK(model->PauliFactor(free.Momentum().P()) == 1);
    CHECK(model->PauliFactor(blocked.Momentum().P()) == 0);

    auto no_pauli = MakeModel("{Distribution: RFG, kF: 221, Binding: 25, PauliBlocking: false}");
    CHECK(no_pauli->PauliFactor(blocked.Momentum().P()) == 1);

    auto cfg = MakeModel(kCFG);
    CHECK_THAT(cfg->PauliFactor(100), WithinRel(1 - cfg->Occupation(100), 1e-12));
    CHECK(cfg->PauliFactor(100) > 0); // partial blocking inside the Fermi sea
}

TEST_CASE("MomentumMapper", "[QEFermiGas]") {
    achilles::ProcessInfo info({achilles::PID::electron(), {achilles::PID::electron()}});
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
    const double mN = achilles::Constant::mN;
    auto model = MakeModel("{Distribution: RFG, kF: 221, Binding: 25}");
    auto mapper = achilles::MomentumMapper::Construct(info, 1);
    mapper->Configure(model->PhaseSpaceSettings());
    mapper->SetMasses({0, mN * mN});
    REQUIRE(mapper->NDims() == 3);

    SECTION("Round trip and on-shell condition") {
        std::mt19937 gen(1234);
        std::uniform_real_distribution<double> uni(0, 1);
        for(const double ebeam : {300.0, 1000.0}) {
            for(size_t i = 0; i < 1000; ++i) {
                std::vector<achilles::FourVector> mom = {{ebeam, 0, 0, ebeam}, {}};
                std::vector<double> ran = {uni(gen), uni(gen), uni(gen)};
                mapper->GeneratePoint(mom, ran);
                std::vector<double> ran2(3);
                const double wgt = mapper->GenerateWeight(mom, ran2);
                for(size_t j = 0; j < 3; ++j) CHECK_THAT(ran2[j], WithinRel(ran[j], 1e-8));
                const double p = mom[1].P();
                CHECK_THAT(mom[1].E(), WithinRel(std::sqrt(p * p + mN * mN) - 25, 1e-12));
                CHECK((mom[0] + mom[1]).M2() >= mN * mN * (1 - 1e-12)); // Smin = mN^2
                CHECK(wgt > 0);
            }
        }
    }

    SECTION("1/density integrates the Fermi sphere volume") {
        // E[1/density] = volume of sampled region = 4/3 pi kF^3 when no kinematic cut applies
        std::mt19937 gen(42);
        std::uniform_real_distribution<double> uni(0, 1);
        constexpr size_t n = 200000;
        double sum = 0;
        for(size_t i = 0; i < n; ++i) {
            std::vector<achilles::FourVector> mom = {{2000, 0, 0, 2000}, {}};
            std::vector<double> ran = {uni(gen), uni(gen), uni(gen)};
            mapper->GeneratePoint(mom, ran);
            std::vector<double> ran2(3);
            sum += 1.0 / mapper->GenerateWeight(mom, ran2);
        }
        CHECK_THAT(sum / n, WithinRel(4.0 / 3 * M_PI * std::pow(221.0, 3), 2e-2));
    }
}
