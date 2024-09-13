#include "Achilles/Particle.hh"
#include "catch2/catch.hpp"

#include "Achilles/CascadeInteractions/DeltaInteractions.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Units.hh"
#include "Achilles/Utilities.hh"

#include <fstream>

using achilles::operator""_GeV;

// TEST_CASE("NNInelasticBkg", "[Delta]") {
//     auto sqrts_vec = achilles::Linspace(2_GeV, 5_GeV, 100);
//     achilles::DeltaInteraction delta;
//     fmt::print("sqrts,xsec_pp,xsec_pn,xsec_pppi0,xsec_pnpip,xsec_pppim,xsec_pnpi0\n");
//     for(const auto &sqrts : sqrts_vec) {
//         auto nn_results = delta.TestNNElastic(sqrts);
//         fmt::print("{:.3e},{:.3e},{:.3e},{:.3e},{:.3e},{:.3e},{:.3e}\n", sqrts, nn_results.first,
//                    nn_results.second, delta.TestInelastic1Pi(sqrts, 0),
//                    delta.TestInelastic1Pi(sqrts, 1), delta.TestInelastic1Pi(sqrts, 2),
//                    delta.TestInelastic1Pi(sqrts, 3));
//     }
// }

TEST_CASE("NN -> NDelta vs. NDelta -> NN", "[Delta]") {
    auto sqrts_vec = achilles::Linspace(1.0, 4.0, 301);
    achilles::DeltaInteraction delta;
    std::ofstream file("delta.txt");
    file << fmt::format("sqrts,dsigmadm,sigmapp,sigmapd\n");

    double expected = 0.2456636303;
    auto result = delta.TestDeltaDSigmaDOmegaDM(-0.81, 2.041, 1.093, achilles::PID::delta0());
    CHECK_THAT(result, Catch::Matchers::WithinAbs(expected, 1e-6));

    for(const auto &sqrts : sqrts_vec) {
        file << fmt::format("{:.3e},{:.3e},{:.3e},{:.3e}\n", sqrts,
                            delta.TestDeltaDSigma(0, sqrts, 1.232), delta.TestDeltaSigma(sqrts),
                            delta.TestDeltaDSigma(1, sqrts, 1.232));
    }

    auto mass_vec = achilles::Linspace(0.0, 4.0, 301);
    std::ofstream file2("delta_mass.txt");
    file2 << fmt::format("mass,dsigmadm\n");

    for(const auto &mass : mass_vec) {
        file2 << fmt::format("{:.3e},{:.3e}\n", mass, delta.TestDeltaDSigma(0, 3.0, mass));
    }
}

TEST_CASE("P Pi+ -> Delta++", "[Delta]") {
    auto plab_vec = achilles::Linspace(0.0, 2.0, 201);
    achilles::DeltaInteraction delta;
    std::ofstream file("ppip.txt");
    file << fmt::format("plab,sigma\n");

    achilles::Particle p1(2212, {achilles::ParticleInfo(2212).Mass(), 0, 0, 0});
    achilles::Particle p2(211);
    achilles::PID delta_id(2224);
    for(const auto &plab : plab_vec) {
        p2.Momentum() = {sqrt(pow(plab * 1_GeV, 2) + pow(p2.Mass(), 2)), 0, 0, plab * 1_GeV};
        file << fmt::format("{:.3e},{:.3e}\n", plab, delta.TestNPiSigma(p1, p2, delta_id));
    }
}

TEST_CASE("P Pi0 -> Delta+", "[Delta]") {
    auto plab_vec = achilles::Linspace(0.0, 2.0, 201);
    achilles::DeltaInteraction delta;
    std::ofstream file("ppi0.txt");
    file << fmt::format("plab,sigma\n");

    achilles::Particle p1(2212, {achilles::ParticleInfo(2212).Mass(), 0, 0, 0});
    achilles::Particle p2(111);
    achilles::PID delta_id(2214);
    for(const auto &plab : plab_vec) {
        p2.Momentum() = {sqrt(pow(plab * 1_GeV, 2) + pow(p2.Mass(), 2)), 0, 0, plab * 1_GeV};
        file << fmt::format("{:.3e},{:.3e}\n", plab, delta.TestNPiSigma(p1, p2, delta_id));
    }
}

TEST_CASE("P Pi- -> Delta0", "[Delta]") {
    auto plab_vec = achilles::Linspace(0.0, 2.0, 201);
    achilles::DeltaInteraction delta;
    std::ofstream file("ppim.txt");
    file << fmt::format("plab,sigma\n");

    achilles::Particle p1(2212, {achilles::ParticleInfo(2212).Mass(), 0, 0, 0});
    achilles::Particle p2(-211);
    achilles::PID delta_id(2114);
    for(const auto &plab : plab_vec) {
        p2.Momentum() = {sqrt(pow(plab * 1_GeV, 2) + pow(p2.Mass(), 2)), 0, 0, plab * 1_GeV};
        file << fmt::format("{:.3e},{:.3e}\n", plab, delta.TestNPiSigma(p1, p2, delta_id));
    }
}

TEST_CASE("TestInterp", "[DeltaSpeed]") {
    achilles::DeltaInteraction delta;
    delta.TestInterpolation();
}
