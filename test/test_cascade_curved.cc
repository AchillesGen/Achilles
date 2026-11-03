// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/generators/catch_generators.hpp"

#include "Achilles/Cascade.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Potential.hh"
#include "Achilles/PotentialFlow.hh"
#include "Achilles/SymplecticIntegrator.hh"

#include <cmath>
#include <memory>

namespace {

/// Central, momentum-independent Woods-Saxon vector potential, V(r) = -V0 / (1 + e^{(r-R)/a})
class WoodsSaxonTestPotential : public achilles::Potential {
  public:
    WoodsSaxonTestPotential(double v0 = 45, double r0 = 2.6, double a = 0.5)
        : m_v0{v0}, m_r0{r0}, m_a{a} {}
    achilles::PotentialVals operator()(double, double r) const override {
        achilles::PotentialVals vals{};
        vals.rvector = -m_v0 / (1 + std::exp((r - m_r0) / m_a));
        return vals;
    }

  private:
    double m_v0, m_r0, m_a;
};

} // namespace

TEST_CASE("PotentialProp force is central", "[CascadeCurved]") {
    auto potential = std::make_shared<WoodsSaxonTestPotential>();
    // Position along x, momentum along y: a central force must point along x, not along y
    const achilles::ThreeVector q{2.0, 0, 0}, p{0, 500, 0};
    const achilles::Particle part{achilles::PID::proton(),
                                  {std::sqrt(500 * 500 + achilles::Constant::mN * achilles::Constant::mN), 0,
                                   500, 0},
                                  q,
                                  achilles::ParticleStatus::propagating};
    auto integrator = achilles::Cascade::MakeIntegrator(part, potential);
    const auto dHdr = integrator.dHdr()(q, p, potential);

    SECTION("dH/dr is along r-hat") {
        CHECK(std::abs(dHdr.Py()) < 1e-10 * dHdr.Magnitude());
        CHECK(std::abs(dHdr.Pz()) < 1e-10 * dHdr.Magnitude());
        // Attractive well rising with r: dH/dr points outward (force inward)
        CHECK(dHdr.Px() > 0);
    }

    SECTION("dH/dr equals the finite-difference gradient of H") {
        constexpr double h = 1e-4;
        const double mass = achilles::Constant::mN;
        const double dHdx = (achilles::flow::Hamiltonian(q + achilles::ThreeVector{h, 0, 0}, p,
                                                         *potential, mass) -
                             achilles::flow::Hamiltonian(q - achilles::ThreeVector{h, 0, 0}, p,
                                                         *potential, mass)) /
                            (2 * h);
        CHECK(dHdr.Px() == Catch::Approx(dHdx).epsilon(1e-6));
    }

    SECTION("Angular momentum is conserved along the integrated path") {
        // Tao's extended phase space conserves L only up to a bounded O(h^2) error: 2.5e-6
        // relative here (2.7e-4 with the original coupling scale), 9.4e-3 with dH/dr along p-hat
        const auto L0 = q.Cross(p);
        for(size_t i = 0; i < 200; ++i) integrator.Step<2>(0.05);
        const auto L1 = integrator.Q().Cross(integrator.P());
        CHECK((L1 - L0).Magnitude() < 1e-3 * L0.Magnitude());
    }
}

TEST_CASE("Tao coupling is dimensionally consistent: steep wall", "[CascadeCurved]") {
    // Sharp well (diffuseness 0.05 fm): with the original coupling (scale 1, fm mixed with MeV)
    // the interior momentum came out 470.5 instead of 531.35 MeV/c and |H - H0| reached 38 MeV
    // at h = 0.04 fm/c, for any step size; with the scale hbar c the flow is reproduced
    auto potential = std::make_shared<WoodsSaxonTestPotential>(40, 6.0, 0.05);
    const double mass = 938.272;
    const achilles::ThreeVector q0{3, 0, -6.8775}, p0{0, 0, 444.6};
    auto dr = [mass](const achilles::ThreeVector &q, const achilles::ThreeVector &p,
                     std::shared_ptr<achilles::Potential> pot) {
        return achilles::flow::dHdr(q, p, *pot, mass);
    };
    auto dp = [mass](const achilles::ThreeVector &q, const achilles::ThreeVector &p,
                     std::shared_ptr<achilles::Potential> pot) {
        return achilles::flow::dHdp(q, p, *pot, mass);
    };
    achilles::SymplecticIntegrator integrator(q0, p0, potential, dr, dp, 20);
    REQUIRE(integrator.CouplingScale() == Catch::Approx(197.3269804));
    const double energy = achilles::flow::Hamiltonian(q0, p0, *potential, mass);
    // inside the flat bottom: E + V0 = sqrt(k^2 + m^2)
    const double kin = std::sqrt(std::pow(energy + 40, 2) - mass * mass);
    double rmin = 1e9, pin = 0, drift = 0;
    for(size_t i = 0; i < 1000; ++i) {
        integrator.Step<2>(0.04);
        const auto &state = integrator.State();
        const auto q = 0.5 * (state.q + state.x), p = 0.5 * (state.p + state.y);
        drift = std::max(drift, std::abs(achilles::flow::Hamiltonian(q, p, *potential, mass) - energy));
        if(q.P() < rmin) {
            rmin = q.P();
            pin = p.P();
        }
    }
    CHECK(pin == Catch::Approx(kin).epsilon(1e-3));
    CHECK(rmin == Catch::Approx(3 * 444.6 / kin).epsilon(1e-3)); // angular momentum b p = rmin k
    CHECK(drift < 0.05);
}

// ---------------------------------------------------------------------------------------------
// Curved-path algorithms
// ---------------------------------------------------------------------------------------------

#include "Achilles/CascadeInteractions/ConstantInteractions.hh"
#include "Achilles/CurvedPath.hh"
#include "Achilles/Event.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Random.hh"
#include "mock_classes.hh"

namespace {

static const double cRadius = 8.0;

class ZeroTestPotential : public achilles::Potential {
  public:
    achilles::PotentialVals operator()(double, double) const override { return {}; }
};

/// Integral of kappa along a straight line at impact parameter b (z = b sinh u)
double LineIntegral(const achilles::PairKernel &kernel, double b, double sigma) {
    constexpr size_t n = 20000;
    constexpr double umax = 12;
    double sum = 0;
    for(size_t i = 0; i <= n; ++i) {
        const double u = umax * static_cast<double>(i) / n;
        const double w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
        const double z = b * std::sinh(u);
        sum += w * kernel.Kappa(std::hypot(b, z), sigma) * b * std::cosh(u);
    }
    return 2 * sum * umax / n / 3;
}

achilles::InteractionHandler ConstantPN(double xsec_mb) {
    achilles::InteractionHandler handler;
    achilles::ConstantInteraction interaction;
    interaction.AddInteraction({achilles::PID::proton(), achilles::PID::neutron()},
                               {{{achilles::PID::proton(), achilles::PID::neutron()}, xsec_mb}});
    interaction.AddInteraction({achilles::PID::neutron(), achilles::PID::neutron()},
                               {{{achilles::PID::neutron(), achilles::PID::neutron()}, xsec_mb}});
    interaction.AddInteraction({achilles::PID::proton(), achilles::PID::proton()},
                               {{{achilles::PID::proton(), achilles::PID::proton()}, xsec_mb}});
    handler.RegisterInteraction(std::make_unique<achilles::ConstantInteraction>(interaction));
    return handler;
}

achilles::Particle Neutron(double x, double y, double z) {
    return achilles::Particle{achilles::PID::neutron(),
                              {achilles::Constant::mN, 0, 0, 0},
                              {x, y, z},
                              achilles::ParticleStatus::background};
}

achilles::Particle Proton(const achilles::ThreeVector &q, const achilles::ThreeVector &p) {
    const double mass = achilles::ParticleInfo(achilles::PID::proton()).Mass();
    return achilles::Particle{achilles::PID::proton(),
                              {std::sqrt(p.P2() + mass * mass), p.Px(), p.Py(), p.Pz()},
                              q,
                              achilles::ParticleStatus::propagating};
}

} // namespace

TEST_CASE("Pair kernel: inverse Abel transform of the Achilles profiles", "[CascadeCurved]") {
    using Profile = achilles::PairKernel::Profile;
    // Reference values from the closed forms (mpmath polylog, scipy k0), sigma = 4 fm^2
    CHECK(achilles::PairKernel(Profile::Gaussian).Kappa(0.5, 4) ==
          Catch::Approx(1.289981062836434).epsilon(1e-10));
    CHECK(achilles::PairKernel(Profile::Gaussian).Kappa(2.0, 4) ==
          Catch::Approx(0.022291402703802816).epsilon(1e-10));
    CHECK(achilles::PairKernel(Profile::Pion).Kappa(0.5, 4) ==
          Catch::Approx(0.5158626415677963).epsilon(1e-8));
    CHECK(achilles::PairKernel(Profile::Pion).Kappa(0.05, 4) ==
          Catch::Approx(9.05773637592298).epsilon(1e-8));

    auto profile = GENERATE(Profile::Gaussian, Profile::Pion);
    const achilles::PairKernel kernel(profile);
    const double sigma = GENERATE(1.0, 4.0);
    for(double b : {0.1, 0.5, 1.2}) {
        // The straight-line integral of kappa reproduces X(b) = -ln(1 - P)
        CHECK(LineIntegral(kernel, b, sigma) ==
              Catch::Approx(kernel.Exponent(b * b, sigma)).epsilon(1e-6));
    }
    // kappa ~ 1/r (Gaussian) and 1/(2r) (Pion) at small r
    const double r = 1e-4;
    CHECK(r * kernel.Kappa(r, sigma) ==
          Catch::Approx(profile == Profile::Gaussian ? 1.0 : 0.5).epsilon(1e-2));
}

TEST_CASE("Curved path: dense output and bounds", "[CascadeCurved]") {
    const double mass = achilles::Constant::mN;
    const achilles::ThreeVector q0{-4, 1.2, 0.3}, p0{340, 20, -10};

    SECTION("V = 0: the integrator gives an exact straight line") {
        auto zero = std::make_shared<ZeroTestPotential>();
        achilles::CurvedPath path(0, q0, p0, mass, zero, 0.2, 20, false);
        const auto v = p0 / std::sqrt(p0.P2() + mass * mass);
        for(double t : {0.05, 1.37, 7.9, 23.3}) {
            CHECK((path.Position(t) - (q0 + v * t)).P() < 1e-12);
            CHECK((path.Velocity(t) - v).P() < 1e-12);
            CHECK((path.Momentum(t) - p0).P() < 1e-9);
        }
    }

    SECTION("Woods-Saxon: bounds hold on the dense path, energy is conserved") {
        auto ws = std::make_shared<WoodsSaxonTestPotential>();
        constexpr double step = 0.04;
        achilles::CurvedPath path(0, q0, p0, mass, ws, step, 8, false);
        const double texit = path.ExitTime();
        REQUIRE(std::isfinite(texit));
        double drift = 0;
        for(size_t k = 0; path.SegmentStart(k) < texit; ++k) {
            const double ta = path.SegmentStart(k);
            const auto bound = path.Bounds(ta, path.SegmentEnd(k));
            for(size_t i = 0; i <= 16; ++i) {
                const double t = ta + step * static_cast<double>(i) / 16;
                // the chord-tube bound: the path stays within delta of the chord
                const auto q = path.Position(t);
                const auto ab = bound.b - bound.a;
                const double s = std::clamp((q - bound.a).Dot(ab) / ab.P2(), 0.0, 1.0);
                CHECK((q - (bound.a + s * ab)).P() <= bound.delta + 1e-12);
                CHECK(path.Velocity(t).P() <= bound.vmax + 1e-12);
            }
            drift = std::max(drift, std::abs(path.Hamiltonian(path.Position(ta), path.Momentum(ta)) -
                                             path.Energy()));
        }
        CHECK(drift < 1e-3); // MeV, at h = 0.04 fm/c (averaged Tao copies)
        // block bounds: every point of a block lies within delta of the block chord
        for(size_t level = 1; level < 6; ++level) {
            for(size_t index = 0; (index + 1) << level <= path.Segment(texit); ++index) {
                const auto &block = path.Block(level, index);
                for(size_t i = 0; i <= 64; ++i) {
                    const double t = block.ta + (block.tb - block.ta) * static_cast<double>(i) / 64;
                    const auto q = path.Position(t);
                    const auto ab = block.b - block.a;
                    const double s = std::clamp((q - block.a).Dot(ab) / ab.P2(), 0.0, 1.0);
                    CHECK((q - (block.a + s * ab)).P() <= block.delta + 1e-12);
                    CHECK(path.Velocity(t).P() <= block.vmax + 1e-12);
                    CHECK(path.Velocity(t).P() >= block.vmin - 1e-12);
                }
            }
        }
        // the path is bent by the attractive well
        CHECK(path.Velocity(texit).Unit().Dot(p0.Unit()) < 0.999);
    }
}

namespace {

struct PairSetup {
    achilles::Cascade::Algorithm alg;
    achilles::Cascade::CurvedOptions::Majorant majorant;
    bool potential;
};

/// Simulate first-collision times of one proton through a few neutrons and compare with
/// SurvivalProfile (no Pauli blocking: hydrogen); also returns the total collision fraction
void CheckSurvival(const PairSetup &setup, const achilles::Particles &targets,
                   const achilles::Particle &projectile, const std::vector<double> &times,
                   size_t ntrials) {
    MockNucleus nucleus;
    auto ws = std::make_shared<WoodsSaxonTestPotential>();
    std::shared_ptr<achilles::Potential> pot = ws;
    ALLOW_CALL(nucleus, Radius()).LR_RETURN((cRadius));
    ALLOW_CALL(nucleus, ID()).RETURN(achilles::PID::hydrogen());
    ALLOW_CALL(nucleus, GetPotential()).LR_RETURN((pot));

    achilles::Cascade cascade(ConstantPN(10), achilles::Cascade::ProbabilityType::Gaussian,
                              setup.alg, achilles::Cascade::InMedium::None,
                              achilles::default_decay_file, setup.potential);
    achilles::Cascade::CurvedOptions opts;
    opts.majorant = setup.majorant;
    cascade.SetCurvedOptions(opts);

    achilles::Event event;
    event.Hadrons() = {projectile};
    event.Hadrons().insert(event.Hadrons().end(), targets.begin(), targets.end());
    const auto expected = cascade.SurvivalProfile(event, &nucleus, 0, times);

    std::vector<double> survived(times.size(), 0);
    for(size_t i = 0; i < ntrials; ++i) {
        event = achilles::Event{};
        event.Hadrons() = {projectile};
        event.Hadrons().insert(event.Hadrons().end(), targets.begin(), targets.end());
        cascade.Evolve(event, &nucleus);
        const auto &diag = cascade.GetDiagnostics();
        CHECK(diag.violations == 0);
        for(size_t k = 0; k < times.size(); ++k)
            survived[k] += !diag.first_set || diag.first_time > times[k];
    }
    CHECK(expected.back() < 0.9);
    for(size_t k = 0; k < times.size(); ++k) {
        const double frac = survived[k] / static_cast<double>(ntrials);
        const double err = std::sqrt(expected[k] * (1 - expected[k]) / static_cast<double>(ntrials));
        INFO("t = " << times[k] << ", simulated " << frac << ", expected " << expected[k]);
        CHECK(frac == Catch::Approx(expected[k]).margin(5 * err + 1e-9));
    }
}

} // namespace

TEST_CASE("Curved algorithms: single pair interacts with P(b^2, sigma) at V = 0",
          "[CascadeCurved]") {
    using Alg = achilles::Cascade::Algorithm;
    const double b = GENERATE(0.3, 0.7);
    auto alg = GENERATE(Alg::Hazard, Alg::CurvedVeto);
    const double sigma = 1; // 10 mb
    const double expected = std::exp(-M_PI * b * b / sigma);
    MockNucleus nucleus;
    ALLOW_CALL(nucleus, Radius()).LR_RETURN((cRadius));
    ALLOW_CALL(nucleus, ID()).RETURN(achilles::PID::hydrogen());
    ALLOW_CALL(nucleus, GetPotential()).RETURN(nullptr);
    achilles::Cascade cascade(ConstantPN(10), achilles::Cascade::ProbabilityType::Gaussian, alg,
                              achilles::Cascade::InMedium::None);
    constexpr size_t ntrials = 20000;
    size_t hits = 0;
    // Start 4 fm before the target so that the whole pair hazard lies ahead
    for(size_t i = 0; i < ntrials; ++i) {
        achilles::Event event;
        event.Hadrons() = {Proton({-4, 0, 0}, {500, 0, 0}), Neutron(0, b, 0)};
        cascade.Evolve(event, &nucleus);
        hits += event.Hadrons()[1].Status() == achilles::ParticleStatus::interacted;
    }
    const double frac = static_cast<double>(hits) / ntrials;
    const double err = std::sqrt(expected * (1 - expected) / ntrials);
    CHECK(frac == Catch::Approx(expected).margin(5 * err));
}

TEST_CASE("Curved algorithms: first-collision times follow SurvivalProfile", "[CascadeCurved]") {
    using Alg = achilles::Cascade::Algorithm;
    using Maj = achilles::Cascade::CurvedOptions::Majorant;
    const achilles::Particles targets = {Neutron(-1.5, 0.9, 0.1), Neutron(0.2, 0.4, -0.2),
                                         Neutron(1.6, -0.3, 0.5)};
    const auto projectile = Proton({-4.5, 1.0, 0}, {340, 0, 0});
    const std::vector<double> times = {4, 8, 12, 16, 20, 30};
    auto setup = GENERATE(PairSetup{Alg::Hazard, Maj::PerStep, false},
                          PairSetup{Alg::Hazard, Maj::PerStep, true},
                          PairSetup{Alg::Hazard, Maj::Constant, true},
                          PairSetup{Alg::Hazard, Maj::Inflated, true},
                          PairSetup{Alg::CurvedVeto, Maj::PerStep, false},
                          PairSetup{Alg::CurvedVeto, Maj::PerStep, true});
    CheckSurvival(setup, targets, projectile, times, 4000);
}

TEST_CASE("CurvedVeto at V = 0 offers the trials of Algorithm: Veto", "[CascadeCurved]") {
    // Independent product over straight-line closest approaches (Particle.cc ClosestApproach)
    MockNucleus nucleus;
    ALLOW_CALL(nucleus, Radius()).LR_RETURN((cRadius));
    ALLOW_CALL(nucleus, ID()).RETURN(achilles::PID::hydrogen());
    ALLOW_CALL(nucleus, GetPotential()).RETURN(nullptr);
    achilles::Cascade cascade(ConstantPN(30), achilles::Cascade::ProbabilityType::Gaussian,
                              achilles::Cascade::Algorithm::CurvedVeto,
                              achilles::Cascade::InMedium::None);
    auto uni = [](double a, double b) { return achilles::Random::Instance().Uniform(a, b); };
    const double sigma = 3; // fm^2
    for(size_t trial = 0; trial < 200; ++trial) {
        achilles::Event event;
        const achilles::ThreeVector p{uni(-400.0, 400.0), uni(-400.0, 400.0),
                                      uni(-400.0, 400.0)};
        event.Hadrons() = {Proton({uni(-2.0, 2.0), uni(-2.0, 2.0), 0}, p)};
        for(size_t j = 0; j < 11; ++j)
            event.Hadrons().push_back(Neutron(uni(-3.0, 3.0), uni(-3.0, 3.0),
                                              uni(-3.0, 3.0)));
        const std::vector<double> times = {0.5, 1, 2, 4, 8, 16, 32};
        const auto profile = cascade.SurvivalProfile(event, &nucleus, 0, times);
        const auto &proj = event.Hadrons()[0];
        const auto v = proj.Beta();
        const double radius2 = cRadius * cRadius;
        for(size_t k = 0; k < times.size(); ++k) {
            double expected = 1;
            for(size_t j = 1; j < event.Hadrons().size(); ++j) {
                const double tc = ClosestApproach(proj, event.Hadrons()[j]);
                const auto at = proj.Position() + v * tc;
                if(tc <= 0 || tc > times[k] || at.P2() > radius2) continue;
                const double b2 = (at - event.Hadrons()[j].Position()).P2();
                expected *= 1 - std::exp(-M_PI * b2 / sigma);
            }
            CHECK(profile[k] == Catch::Approx(expected).margin(1e-12));
        }
    }
}
