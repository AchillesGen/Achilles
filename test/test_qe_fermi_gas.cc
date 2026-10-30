// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/Nucleus.hh"
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

namespace {

// Radial bins of width dr out to rmax with densities rho_p(r), rho_n(r) at the bin centres
template <class Fp, class Fn>
achilles::LocalFermiSea MakeSea(double rmax, Fp rho_p, Fn rho_n, double dr = 0.005) {
    const auto n = static_cast<size_t>(std::round(rmax / dr));
    std::vector<double> edges(n + 1), rp(n), rn(n);
    for(size_t i = 0; i <= n; ++i) edges[i] = static_cast<double>(i) * dr;
    for(size_t i = 0; i < n; ++i) {
        const double r = 0.5 * (edges[i] + edges[i + 1]);
        rp[i] = rho_p(r);
        rn[i] = rho_n(r);
    }
    return {edges, rp, rn};
}

// 4 pi int p^2 dp 2 V(p) / (2 pi hbar c)^3: the number of nucleons in the local Fermi sea
double LocalNormalization(const achilles::LocalFermiSea &sea, achilles::PID nucleon) {
    constexpr size_t n = 20000;
    const double h = 1.01 * sea.MaxFermiMomentum() / n;
    double sum = 0;
    for(size_t i = 0; i < n; ++i) {
        const double p = (static_cast<double>(i) + 0.5) * h;
        sum += 4 * M_PI * p * p * 2 * sea.Volume(p, nucleon) * h;
    }
    return sum / std::pow(2 * M_PI * achilles::Constant::HBARC, 3);
}

} // namespace

TEST_CASE("Local Fermi sea", "[QEFermiGas]") {
    const auto proton = achilles::PID::proton(), neutron = achilles::PID::neutron();

    SECTION("A uniform sphere reproduces the global Fermi gas") {
        constexpr double R0 = 2.5, Z = 6;
        const double rho0 = Z / (4.0 / 3 * M_PI * std::pow(R0, 3));
        const double kf = achilles::Constant::HBARC * std::cbrt(3 * M_PI * M_PI * rho0);
        auto uniform = [&](double r) { return r < R0 ? rho0 : 0.0; };
        auto sea = MakeSea(R0, uniform, uniform);
        const double v0 = 4.0 / 3 * M_PI * std::pow(R0, 3);

        CHECK_THAT(sea.MaxFermiMomentum(), WithinRel(kf, 1e-12));
        CHECK_THAT(sea.Volume(0.5 * kf, proton), WithinRel(v0, 1e-10));
        CHECK(sea.Volume(1.01 * kf, proton) == 0);
        // 2 V0 / (2 pi hbar c)^3 = Z * 3 / (4 pi kF^3): the RFG weight per d^3p
        CHECK_THAT(2 * v0 / std::pow(2 * M_PI * achilles::Constant::HBARC, 3),
                   WithinRel(Z * 3 / (4 * M_PI * std::pow(kf, 3)), 1e-12));
        CHECK(sea.PauliFactor(0.5 * kf, proton, 1.1 * kf, proton) == 1);
        CHECK(sea.PauliFactor(0.5 * kf, proton, 0.9 * kf, proton) == 0);
    }

    // Carbon-like profile with a central depression, as in the 12C VMC density
    auto dip = [](double r) {
        return 0.09 * (1 - 0.4 * std::exp(-r * r / 0.8)) / (1 + std::exp((r - 2.3) / 0.45));
    };
    auto sea = MakeSea(8, dip, [&](double r) { return 1.2 * dip(r); });

    SECTION("The momentum integral counts every nucleon") {
        CHECK_THAT(LocalNormalization(sea, proton), WithinRel(sea.NNucleons(proton), 2e-3));
        CHECK_THAT(LocalNormalization(sea, neutron), WithinRel(sea.NNucleons(neutron), 2e-3));
    }

    SECTION("V(k) falls with k and the Pauli factor is a fraction of it") {
        double last = sea.Volume(0, proton);
        for(double k = 10; k < sea.MaxFermiMomentum(); k += 10) {
            const double v = sea.Volume(k, proton);
            CHECK(v <= last);
            last = v;
        }
        for(double p = 20; p < 260; p += 40) {
            for(double pout = 0; pout < 400; pout += 25) {
                const double f = sea.PauliFactor(p, proton, pout, proton);
                CHECK(f >= 0);
                CHECK(f <= 1);
                if(pout <= p) CHECK(f == 0);
                // Same nucleon type: 1 - V(p') / V(p), also with a non-monotonic density
                if(pout > p && sea.Volume(p, proton) > 0)
                    CHECK_THAT(
                        f, WithinRel(1 - sea.Volume(pout, proton) / sea.Volume(p, proton), 1e-12) ||
                               Catch::Matchers::WithinAbs(0, 1e-12));
            }
        }
    }

    SECTION("The allowed radii enclose every bin where the struck nucleon can sit") {
        constexpr double dr = 0.005; // bin width used by MakeSea
        auto allowed = [&](double r, double p, double pout) {
            const auto kf = [&](double scale) {
                return achilles::Constant::HBARC * std::cbrt(3 * M_PI * M_PI * scale * dip(r));
            };
            return kf(1.2) > p && kf(1.0) <= pout; // neutron occupied, proton state free
        };
        for(double p = 20; p < 260; p += 40) {
            for(double pout = p + 10; pout < 300; pout += 30) {
                const auto radii = sea.AllowedRadii(p, neutron, pout, proton);
                if(sea.PauliFactor(p, neutron, pout, proton) <= 0) {
                    CHECK_FALSE(radii);
                    continue;
                }
                REQUIRE(radii);
                CHECK(allowed(radii->first + dr / 2, p, pout));
                CHECK(allowed(radii->second - dr / 2, p, pout));
                for(double r = dr / 2; r < 8; r += dr)
                    if(r < radii->first || r > radii->second) CHECK_FALSE(allowed(r, p, pout));
            }
        }
        CHECK_FALSE(sea.AllowedRadii(100, proton, 90, proton));

        // Off-centre peak: a momentum above the central kF is only allowed away from the centre
        const double kf0 = achilles::Constant::HBARC * std::cbrt(3 * M_PI * M_PI * dip(0));
        const auto off_centre = sea.AllowedRadii(kf0 + 10, proton, 400, proton);
        REQUIRE(off_centre);
        CHECK(off_centre->first > 0);

        constexpr double R0 = 2.5;
        auto uniform = [&](double r) { return r < R0 ? 0.1 : 0.0; };
        const auto ball = MakeSea(R0, uniform, uniform).AllowedRadii(50, proton, 400, proton);
        REQUIRE(ball);
        CHECK(ball->first == 0);
        CHECK_THAT(ball->second, WithinRel(R0, 1e-12));
    }

    SECTION("A neutron turning into a proton is blocked by the proton sea") {
        // Neutron sea deeper than the proton sea: p' between the two kF is free
        const double kfp = sea.MaxFermiMomentum() / std::cbrt(1.2);
        CHECK(sea.PauliFactor(10, neutron, 1.01 * kfp, proton) == 1);
        CHECK(sea.PauliFactor(10, neutron, 1.01 * kfp, neutron) < 1);
    }
}

TEST_CASE("QEFermiGas with a local Fermi gas", "[QEFermiGas]") {
    SECTION("kF: Local is read from the run card, and only for RFG") {
        auto model = MakeModel("{Distribution: RFG, kF: Local, Binding: 25}");
        CHECK(model->IsLocal());
        CHECK_FALSE(MakeModel("{Distribution: RFG, kF: 221, Binding: 25}")->IsLocal());
        CHECK_THROWS(
            MakeModel("{Distribution: CFG, kF: Local, Binding: 25, c0: 4.16, lambda: 2.75}"));
    }

    SECTION("SetNucleus builds the sea from the nucleus densities") {
        static constexpr size_t Z = 6;
        achilles::Particles particles;
        for(size_t i = 0; i < Z; ++i) {
            particles.emplace_back(achilles::PID::proton());
            particles.emplace_back(achilles::PID::neutron());
        }
        achilles::Nucleus::FermiGas fermiGas = {achilles::Nucleus::FermiGasType::Local, {}};
        auto density = std::make_unique<MockDensity>();
        REQUIRE_CALL(*density, GetConfiguration()).TIMES(1).RETURN(particles);
        achilles::Nucleus nuc(Z, 2 * Z, 0, 225, "data/densities/c12.prova.txt",
                              "data/densities/c12.prova.txt", fermiGas, std::move(density));

        auto model = MakeModel("{Distribution: RFG, kF: Local, Binding: 25}");
        model->SetNucleus(nuc);
        CHECK_THAT(model->Sea().NNucleons(achilles::PID::proton()), WithinRel(6.0, 5e-3));
        // The 12C density peaks off-centre: the largest local kF is not the central one
        CHECK(model->PCut() > nuc.FermiMomentum(0, achilles::PID::proton()));
        CHECK_THAT(model->PhaseSpaceSettings()["PCut"].as<double>(),
                   WithinRel(model->PCut(), 1e-12));

        // Only the local Fermi gas restricts the struck nucleon's radius
        const auto proton = achilles::PID::proton();
        const auto radii = model->StruckNucleonRadii(100, proton, 400, proton);
        REQUIRE(radii);
        CHECK(radii->first < radii->second);
        CHECK(radii->second <= nuc.Radius());
        CHECK_FALSE(MakeModel("{Distribution: RFG, kF: 221, Binding: 25}")
                        ->StruckNucleonRadii(100, proton, 400, proton));
    }
}
