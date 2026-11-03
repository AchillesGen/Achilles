// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Fast unit tests of the Sudakov veto cascade. The high-statistics validation (Step 4 a-e)
// is scripts/sudakov_validation/run_all.sh + analyse.py.

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Sudakov.hh"

#include <cmath>
#include <random>

#include "yaml-cpp/yaml.h"

using namespace achilles;
using namespace achilles::sudakov;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
std::mt19937_64 rng(20261026);
double U() {
    return (static_cast<double>(rng() >> 11) + 0.5) * 0x1.0p-53;
}

constexpr double kR = 3.0, kRho0 = 12.0 / (4.0 / 3.0 * M_PI * 27.0), kSigma = 40.0;

double ClosedT(double mu) {
    const double y = 2 * mu * kR;
    return 3 / (4 * kR * kR * kR) *
           (kR * kR * (1 - std::exp(-y)) / mu -
            (2 - std::exp(-y) * (y * y + 2 * y + 2)) / (4 * mu * mu * mu));
}

Particle Proton(const ThreeVector &x, ParticleStatus st) {
    const double m = ParticleInfo(PID::proton()).Mass(), p = 1000;
    return Particle(PID::proton(), FourVector{std::sqrt(m * m + p * p), 0, 0, p}, x, st);
}

std::shared_ptr<Engine> TrackA(double headroom, std::vector<double> vars = {}) {
    auto coll = std::make_shared<ForwardCollision>(kSigma);
    DensityRateSettings s;
    s.headroom = headroom;
    s.sigma_scale_vars = vars;
    auto occ = std::make_shared<Occupancy>(Occupancy::Shape::Step, 0.05);
    occ->SetDensities([](double) { return 0.0; }, [](double) { return 0.0; });
    auto rate = std::make_shared<DensityRate>(s, coll, occ);
    rate->SetDensities([](double r) { return r < kR ? kRho0 / 2 : 0.0; },
                       [](double r) { return r < kR ? kRho0 / 2 : 0.0; }, 12, kR);
    EngineSettings es;
    for(size_t i = 0; i < vars.size(); ++i) es.variations.push_back("v" + std::to_string(i));
    return std::make_shared<Engine>(rate, coll, es, U);
}

ThreeVector UniformInSphere() {
    ThreeVector x;
    do {
        x = {kR * (2 * U() - 1), kR * (2 * U() - 1), kR * (2 * U() - 1)};
    } while(x.Magnitude2() >= kR * kR);
    return x;
}
} // namespace

TEST_CASE("Sudakov: Ein and majorant inversion", "[Sudakov]") {
    CHECK_THAT(Ein(1.0), WithinAbs(0.7965995992970534, 1e-13));
    CHECK_THAT(Ein(3.0),
               WithinAbs(std::expint(-3.0) * -1.0 + std::log(3.0) + 0.5772156649015329, 1e-13));
    Majorant M;
    M.Add(0, 1, 2.0);
    M.Add(1, 3, 0.5);
    M.Add(4, 5, 1.0);
    // cumulative: 2 at t=1, 3 at t=3, 4 at t=5
    CHECK_THAT(M.NextContinuous(0, std::exp(-1.0)), WithinAbs(0.5, 1e-14));
    CHECK_THAT(M.NextContinuous(0, std::exp(-2.5)), WithinAbs(2.0, 1e-14));
    CHECK_THAT(M.NextContinuous(0, std::exp(-3.5)), WithinAbs(4.5, 1e-14));
    CHECK(std::isinf(M.NextContinuous(0, std::exp(-4.5))));
    CHECK_THAT(M.NextContinuous(2.0, std::exp(-1.0)), WithinAbs(4.5, 1e-14));
    CHECK_THAT(M.Integral(0.5, 4.5), WithinAbs(1.0 + 1.0 + 0.5, 1e-14));
    CHECK_THROWS(M.Add(2, 6, 1.0));
}

TEST_CASE("Sudakov: smooth occupancy is normalised", "[Sudakov]") {
    Occupancy occ(Occupancy::Shape::FermiDirac, 0.08);
    occ.SetDensities([](double) { return 0.08; }, [](double) { return 0.08; });
    // 3 int Phi(x) x^2 dx = 1
    double sum = 0;
    const size_t n = 20000;
    for(size_t i = 0; i < n; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * occ.XMax() / static_cast<double>(n);
        sum += 3 * x * x * occ.Phi(x) * occ.XMax() / static_cast<double>(n);
    }
    CHECK_THAT(sum, WithinAbs(1.0, 1e-4));
    // <k^2> of the sampled targets equals the occupancy moment
    const double kf = occ.KFermi(1.0, PID::proton());
    double m2 = 0, m2ref = 0;
    for(size_t i = 0; i < 40000; ++i) {
        const double k = occ.SampleK(1.0, PID::proton(), U());
        m2 += k * k / 40000;
    }
    for(size_t i = 0; i < n; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * occ.XMax() / static_cast<double>(n);
        m2ref += 3 * x * x * occ.Phi(x) * x * x * kf * kf * occ.XMax() / static_cast<double>(n);
    }
    CHECK_THAT(m2, WithinRel(m2ref, 0.02));
}

TEST_CASE("Sudakov Track A: uniform-sphere transparency (closed form)", "[Sudakov]") {
    auto eng = TrackA(1.5);
    const size_t N = 20000;
    size_t n0 = 0;
    for(size_t ev = 0; ev < N; ++ev) {
        Event event;
        event.Hadrons().push_back(Proton(UniformInSphere(), ParticleStatus::internal_test));
        eng->Evolve(event, nullptr);
        n0 += eng->GetContext().collisions.empty();
    }
    const double T = static_cast<double>(n0) / N, Tex = ClosedT(kSigma * 0.1 * kRho0);
    CHECK(std::abs(T - Tex) < 4 * std::sqrt(Tex * (1 - Tex) / N));
    CHECK(eng->GetDiagnostics().violations == 0);
}

TEST_CASE("Sudakov Track A: reweighting closure and violation counter", "[Sudakov]") {
    auto eng = TrackA(2.0, {0.8});
    const size_t N = 20000;
    double sw = 0, sw2 = 0;
    for(size_t ev = 0; ev < N; ++ev) {
        Event event;
        event.Hadrons().push_back(Proton(UniformInSphere(), ParticleStatus::internal_test));
        eng->Evolve(event, nullptr);
        const double x = eng->GetContext().collisions.empty() ? eng->GetContext().weights[0] : 0.0;
        sw += x;
        sw2 += x * x;
    }
    const double mean = sw / N, err = std::sqrt((sw2 / N - mean * mean) / N);
    CHECK(std::abs(mean - ClosedT(0.8 * kSigma * 0.1 * kRho0)) < 4 * err);

    auto bad = TrackA(0.7);
    for(size_t ev = 0; ev < 200; ++ev) {
        Event event;
        event.Hadrons().push_back(Proton(UniformInSphere(), ParticleStatus::internal_test));
        bad->Evolve(event, nullptr);
    }
    CHECK(bad->GetDiagnostics().violations > 0);
}

TEST_CASE("Sudakov Track B: toy nucleus, exact first-collision law (atom profile)", "[Sudakov]") {
    auto coll = std::make_shared<ForwardCollision>(kSigma);
    ConfigurationRateSettings s;
    s.law = ConfigurationRateSettings::Law::Kappa1;
    s.profile = ConfigurationRateSettings::Profile::Atom;
    auto rate = std::make_shared<ConfigurationRate>(s, coll);
    Engine eng(rate, coll, EngineSettings{}, U);
    const std::vector<ThreeVector> pts{{0.3, 0, -0.4}, {0.5, 0.2, 0.4}};
    const double p1 = rate->PairProb(0.09, kSigma * 0.1), p2 = rate->PairProb(0.29, kSigma * 0.1);
    const size_t N = 50000;
    size_t first1 = 0, none = 0;
    for(size_t ev = 0; ev < N; ++ev) {
        Event event;
        event.Hadrons().push_back(Proton({0, 0, -10}, ParticleStatus::external_test));
        for(const auto &x : pts)
            event.Hadrons().emplace_back(PID::neutron(),
                                         FourVector{ParticleInfo(PID::neutron()).Mass(), 0, 0, 0},
                                         x, ParticleStatus::background);
        eng.Evolve(event, nullptr);
        const auto &c = eng.GetContext().collisions;
        if(c.empty())
            none++;
        else if(c.front().target == 1)
            first1++;
    }
    const double f1 = static_cast<double>(first1) / N, f0 = static_cast<double>(none) / N;
    CHECK(std::abs(f1 - p1) < 4 * std::sqrt(p1 * (1 - p1) / N));
    const double pn = (1 - p1) * (1 - p2);
    CHECK(std::abs(f0 - pn) < 4 * std::sqrt(pn * (1 - pn) / N) + 1e-4);
}

TEST_CASE("Sudakov: coherent-veto form factor and factory options", "[Sudakov]") {
    // uniform sphere of radius R: F(q) = 3 j1(qR) / (qR)
    AchillesCollisionSettings cs;
    cs.coherent = true;
    AchillesCollision coll(std::make_shared<Occupancy>(Occupancy::Shape::Step, 0.05), cs);
    const double R = 2.9;
    auto rho = [R](double r) { return r < R ? 0.08 : 0.0; };
    coll.SetFormFactors(rho, rho);
    CHECK_THAT(coll.FormFactor(PID::proton(), 0), WithinAbs(1.0, 1e-12));
    for(double q : {50.0, 100.0, 150.0, 213.0, 300.0}) {
        const double x = q / Constant::HBARC * R;
        const double exact = 3 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
        CHECK_THAT(coll.FormFactor(PID::neutron(), q), WithinAbs(exact, 2e-3));
    }
    // factory: defaults and validation
    YAML::Node node = YAML::Load("{Pauli: None, CoherentVeto: true}");
    auto eng = MakeEngine(node, 'B', U);
    CHECK(eng->Settings().single_attempt);
    auto *ac = dynamic_cast<AchillesCollision *>(&eng->Coll());
    REQUIRE(ac != nullptr);
    CHECK_FALSE(ac->Settings().pauli);
    CHECK(ac->Settings().coherent);
    CHECK_FALSE(MakeEngine(YAML::Load("{PauliAttempt: Retry}"), 'A', U)->Settings().single_attempt);
    CHECK_THROWS(MakeEngine(YAML::Load("{Pauli: Off}"), 'A', U));
    CHECK_THROWS(MakeEngine(YAML::Load("{PauliAttempt: Twice}"), 'A', U));
}

TEST_CASE("Sudakov: tabulated coherent veto", "[Sudakov]") {
    auto eng = MakeEngine(YAML::Load("{CoherentTable: {q: [0, 100, 200], v: [1.0, 0.5, 0.0]}}"), 'B', U);
    auto *ac = dynamic_cast<AchillesCollision *>(&eng->Coll());
    REQUIRE(ac != nullptr);
    CHECK(ac->Settings().coherent);
    CHECK_THAT(ac->CoherentVetoProb(PID::proton(), 50), WithinAbs(0.75, 1e-12));
    CHECK_THAT(ac->CoherentVetoProb(PID::neutron(), 150), WithinAbs(0.25, 1e-12));
    CHECK_THAT(ac->CoherentVetoProb(PID::proton(), 900), WithinAbs(0.0, 1e-12));
    CHECK_THROWS(MakeEngine(YAML::Load("{CoherentTable: {q: [0, 100], v: [1.0]}}"), 'B', U));
}
