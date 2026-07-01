#include "Achilles/Particle.hh"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include "Achilles/CascadeInteractions/DeltaInteractions.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Units.hh"
#include "Achilles/Utilities.hh"

#include <algorithm>
#include <cmath>
#include <vector>

using achilles::operator""_GeV;

namespace {

// These tests used to dump cross-section sweeps to text files for visual
// comparison against external calculations. Instead of populating the file
// system, we keep sweeping the same kinematics but assert the results are
// physical (finite and non-negative) and, for resonant channels, that the
// cross section peaks near the expected resonance location.

bool AllPhysical(const std::vector<double> &xsecs) {
    return std::all_of(xsecs.begin(), xsecs.end(),
                       [](double x) { return std::isfinite(x) && x >= 0; });
}

// Index of the largest entry; used to locate a resonance peak in a sweep.
size_t ArgMax(const std::vector<double> &xsecs) {
    return static_cast<size_t>(
        std::distance(xsecs.begin(), std::max_element(xsecs.begin(), xsecs.end())));
}

} // namespace

TEST_CASE("NN -> NDelta vs. NDelta -> NN", "[Delta]") {
    achilles::DeltaInteraction delta;

    double expected = 0.1946713104;
    auto result = delta.TestDeltaDSigmaDOmegaDM(-0.81, 2.041, 1.093, achilles::PID::delta0());
    CHECK_THAT(result, Catch::Matchers::WithinAbs(expected, 1e-6));

    auto sqrts_vec = achilles::Linspace(1.0, 4.0, 301);
    std::vector<double> dsigmadm, sigmapp, sigmapd;
    for(const auto &sqrts : sqrts_vec) {
        dsigmadm.push_back(delta.TestDeltaDSigma(0, sqrts, 1.232));
        sigmapp.push_back(delta.TestDeltaSigma(sqrts));
        sigmapd.push_back(delta.TestDeltaDSigma(1, sqrts, 1.232));
    }
    CHECK(AllPhysical(dsigmadm));
    CHECK(AllPhysical(sigmapp));
    CHECK(AllPhysical(sigmapd));

    auto mass_vec = achilles::Linspace(0.01, 4.0, 301);
    std::vector<double> dsigmadm_vs_mass;
    for(const auto &mass : mass_vec) {
        dsigmadm_vs_mass.push_back(delta.TestDeltaDSigma(0, 3.0, mass));
    }
    CHECK(AllPhysical(dsigmadm_vs_mass));
}

TEST_CASE("P Pi+ -> Delta++", "[Delta]") {
    auto plab_vec = achilles::Linspace(0.0, 2.0, 201);
    achilles::DeltaInteraction delta;

    achilles::Particle p1(2212, {achilles::ParticleInfo(2212).Mass(), 0, 0, 0});
    achilles::Particle p2(211);
    achilles::PID delta_id(2224);
    std::vector<double> sigma;
    for(const auto &plab : plab_vec) {
        p2.Momentum() = {sqrt(pow(plab * 1_GeV, 2) + pow(p2.Mass(), 2)), 0, 0, plab * 1_GeV};
        sigma.push_back(delta.TestNPiSigma(p1, p2, delta_id));
    }

    CHECK(AllPhysical(sigma));
    // The Delta(1232) resonance peaks around plab ~ 0.30 GeV for pi N scattering.
    CHECK_THAT(plab_vec[ArgMax(sigma)], Catch::Matchers::WithinAbs(0.30, 0.1));
}

TEST_CASE("P Pi0 -> Delta+", "[Delta]") {
    auto plab_vec = achilles::Linspace(0.0, 2.0, 201);
    achilles::DeltaInteraction delta;

    achilles::Particle p1(2212, {achilles::ParticleInfo(2212).Mass(), 0, 0, 0});
    achilles::Particle p2(111);
    achilles::PID delta_id(2214);
    std::vector<double> sigma;
    for(const auto &plab : plab_vec) {
        p2.Momentum() = {sqrt(pow(plab * 1_GeV, 2) + pow(p2.Mass(), 2)), 0, 0, plab * 1_GeV};
        sigma.push_back(delta.TestNPiSigma(p1, p2, delta_id));
    }

    CHECK(AllPhysical(sigma));
    CHECK_THAT(plab_vec[ArgMax(sigma)], Catch::Matchers::WithinAbs(0.30, 0.1));
}

TEST_CASE("P Pi- -> Delta0", "[Delta]") {
    auto plab_vec = achilles::Linspace(0.0, 2.0, 201);
    achilles::DeltaInteraction delta;

    achilles::Particle p1(2212, {achilles::ParticleInfo(2212).Mass(), 0, 0, 0});
    achilles::Particle p2(-211);
    achilles::PID delta_id(2114);
    std::vector<double> sigma;
    for(const auto &plab : plab_vec) {
        p2.Momentum() = {sqrt(pow(plab * 1_GeV, 2) + pow(p2.Mass(), 2)), 0, 0, plab * 1_GeV};
        sigma.push_back(delta.TestNPiSigma(p1, p2, delta_id));
    }

    CHECK(AllPhysical(sigma));
    CHECK_THAT(plab_vec[ArgMax(sigma)], Catch::Matchers::WithinAbs(0.30, 0.1));
}

TEST_CASE("P D+ -> N D", "[Delta2]") {
    auto plab_vec = achilles::Linspace(0.01, 2.0, 201);
    achilles::DeltaInteraction delta;

    achilles::Particle p1(2212, {achilles::ParticleInfo(2212).Mass(), 0, 0, 0});
    achilles::Particle p2(2214);
    achilles::PID delta_id1(2214), delta_id2(2224);
    achilles::PID nucleon_id1(2212), nucleon_id2(2112);
    std::vector<double> sigma1, sigma2;
    for(const auto &plab : plab_vec) {
        p2.Momentum() = {sqrt(pow(plab * 1_GeV, 2) + pow(p2.Mass(), 2)), 0, 0, plab * 1_GeV};
        sigma1.push_back(delta.TestNDelta2NDelta(p2, p1, delta_id1, nucleon_id1));
        sigma2.push_back(delta.TestNDelta2NDelta(p2, p1, delta_id2, nucleon_id2));
    }
    CHECK(AllPhysical(sigma1));
    CHECK(AllPhysical(sigma2));

    auto mass_vec = achilles::Linspace(0.938, 2.2, 201);
    double plab = 1.0;
    std::vector<double> sigma1_vs_mass, sigma2_vs_mass;
    for(const auto &mass : mass_vec) {
        p2.Momentum() = {sqrt(pow(plab * 1_GeV, 2) + pow(mass * 1_GeV, 2)), 0, 0, plab * 1_GeV};
        sigma1_vs_mass.push_back(delta.TestNDelta2NDelta(p2, p1, delta_id1, nucleon_id1));
        sigma2_vs_mass.push_back(delta.TestNDelta2NDelta(p2, p1, delta_id2, nucleon_id2));
    }
    CHECK(AllPhysical(sigma1_vs_mass));
    CHECK(AllPhysical(sigma2_vs_mass));
}
