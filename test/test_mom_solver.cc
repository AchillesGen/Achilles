// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/generators/catch_generators_adapters.hpp"
#include "catch2/generators/catch_generators_random.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "mock_classes.hh"

#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/MomSolver.hh"
#include "Achilles/Potential.hh"
#include "Achilles/ThreeVector.hh"
#include "Approx.hh"

#include "catch_utils.hh"
#include "spdlog/spdlog.h"

namespace units = achilles::units;
using namespace achilles::units::literals;

TEST_CASE("CM Momentum Delta solver", "[DeltaFunctions]") {
    achilles::FourVector p1 = achilles::FourVector::FromNative(
        {sqrt(pow(achilles::Constant::mN.native(), 2) + 100 * 100), 0, 0, 100});
    achilles::FourVector p2 = achilles::FourVector::FromNative(
        {sqrt(pow(achilles::Constant::mN.native(), 2) + 100 * 100), 0, 0, -100});
    double cosTheta = 0.5;
    double phi = 0;
    auto p3 = achilles::SolveDelta(p1, p2, achilles::Constant::mN.native(),
                                   achilles::Constant::mN.native(), cosTheta, phi);
    auto p4 = p1 + p2 - p3;

    double sinTheta = sqrt(0.75);
    CHECK(p1 + p2 == p3 + p4);
    CHECK_THAT(p3, FourVectorApprox(achilles::FourVector::FromNative(
                                        {p1.E().native(), 100 * sinTheta, 0, 100 * cosTheta}))
                       .margin(1e-8));
}

TEST_CASE("Arbitrary Frame Delta solver", "[DeltaFunctions]") {
    std::mt19937 m_rand{std::random_device{}()};
    std::uniform_real_distribution<double> m_dist;

    auto p1 = GENERATE(take(10, randomMomentum(1000, achilles::Constant::mN.native())));
    auto p2 = GENERATE(take(10, randomMomentum(1000, achilles::Constant::mN.native())));
    double cosTheta = 2 * m_dist(m_rand) - 1;
    double phi = 2 * M_PI * m_dist(m_rand);

    auto p3 = achilles::SolveDelta(p1, p2, achilles::Constant::mN.native(),
                                   achilles::Constant::mN.native(), cosTheta, phi);
    auto p4 = p1 + p2 - p3;
    auto pcm1 = p1 + p2;
    auto pcm2 = p3 + p4;
    CHECK((pcm1.E().native() == Catch::Approx(pcm2.E().native()) &&
           pcm1.Px().native() == Catch::Approx(pcm2.Px().native()) &&
           pcm1.Py().native() == Catch::Approx(pcm2.Py().native()) &&
           pcm1.Pz().native() == Catch::Approx(pcm2.Pz().native())));
    CHECK(p3.M().native() == Catch::Approx(achilles::Constant::mN.native()));
    CHECK(p4.M().native() == Catch::Approx(achilles::Constant::mN.native()));
}

TEST_CASE("Potential Delta solver", "[DeltaFunctions]") {
    std::mt19937 m_rand{std::random_device{}()};
    std::uniform_real_distribution<double> m_dist;
    constexpr double rmin = 0.0, rmax = 6.0;
    auto p1 = GENERATE(take(10, randomMomentum(1000, achilles::Constant::mN.native())));
    auto p2 = GENERATE(take(10, randomMomentum(225, achilles::Constant::mN.native())));
    auto r1 = GENERATE(take(3, random(rmin, rmax)));
    auto r2 = GENERATE(take(3, random(rmin, rmax)));
    auto q_free = p1 + p2;
    auto nucleus = std::make_shared<MockNucleus>();

    SECTION("WiringaPotential") {
        constexpr double rho0 = 0.16;
        REQUIRE_CALL(*nucleus, Rho(trompeloeil::gt(units::Length{})))
            .LR_RETURN((rho0))
            .TIMES(AT_LEAST(4));

        achilles::WiringaPotential potential(nucleus, rho0);

        auto potential1 = potential(p1.P().native(), r1);
        auto potential2 = potential(p2.P().native(), r2);
        p1.E() =
            units::Energy{sqrt(p1.P2().native() + pow(p1.M().native() + potential1.rscalar, 2)) +
                          potential1.rvector};
        p2.E() =
            units::Energy{sqrt(p2.P2().native() + pow(p2.M().native() + potential2.rscalar, 2)) +
                          potential2.rvector};
        auto q = p1 + p2;

        // Rotate so (p1+p2) to be along z-axis
        auto rotation = q.AlignZ();
        q = q.Rotate(rotation);

        // Randomly generate p3Mag and phi
        double phi = 2 * M_PI * m_dist(m_rand);
        auto range = achilles::FindMomentumRange(q_free);
        const double dp3 = range.second - range.first;
        achilles::FourVector p3;
        while(true) {
            try {
                double p3Mag = dp3 * m_dist(m_rand) + range.first;

                // Solve for p3 and rotate back and calculate p4
                p3 = achilles::SolveDeltaWithPotential(
                    q, potential, achilles::Constant::mN.native(), achilles::Constant::mN.native(),
                    p3Mag, phi, r1, r2);
                break;
            } catch(const std::domain_error &e) { continue; }
        }
        p3 = p3.RotateBack(rotation);
        auto p4 = p1 + p2 - p3;

        auto pcm1 = p1 + p2;
        auto pcm2 = p3 + p4;
        CHECK(pcm1.E().native() == Catch::Approx(pcm2.E().native()));
        CHECK(pcm1.Px().native() == Catch::Approx(pcm2.Px().native()));
        CHECK(pcm1.Py().native() == Catch::Approx(pcm2.Py().native()));
        CHECK(pcm1.Pz().native() == Catch::Approx(pcm2.Pz().native()));
    }

    SECTION("CooperPotential") {
        constexpr size_t AA = 12;
        REQUIRE_CALL(*nucleus, NNucleons()).LR_RETURN((AA)).TIMES(AT_LEAST(4));

        achilles::CooperPotential potential(nucleus);

        auto potential1 = potential(p1.P().native(), 1);
        auto potential2 = potential(p2.P().native(), 1);
        p1.E() =
            units::Energy{sqrt(p1.P2().native() + pow(p1.M().native() + potential1.rscalar, 2)) +
                          potential1.rvector};
        p2.E() =
            units::Energy{sqrt(p2.P2().native() + pow(p2.M().native() + potential2.rscalar, 2)) +
                          potential2.rvector};
        auto q = p1 + p2;

        // Rotate so (p1+p2) to be along z-axis
        auto rotation = q.AlignZ();
        q = q.Rotate(rotation);

        // Randomly generate p3Mag and phi
        double phi = 2 * M_PI * m_dist(m_rand);
        auto range = achilles::FindMomentumRange(q_free);
        const double dp3 = range.second - range.first;
        achilles::FourVector p3;
        while(true) {
            try {
                double p3Mag = dp3 * m_dist(m_rand) + range.first;

                // Solve for p3 and rotate back and calculate p4
                p3 = achilles::SolveDeltaWithPotential(
                    q, potential, achilles::Constant::mN.native(), achilles::Constant::mN.native(),
                    p3Mag, phi, 1, 1);
                break;
            } catch(const std::domain_error &e) { continue; }
        }
        p3 = p3.RotateBack(rotation);
        auto p4 = p1 + p2 - p3;

        auto pcm1 = p1 + p2;
        auto pcm2 = p3 + p4;
        CHECK(pcm1.E().native() == Catch::Approx(pcm2.E().native()));
        CHECK(pcm1.Px().native() == Catch::Approx(pcm2.Px().native()));
        CHECK(pcm1.Py().native() == Catch::Approx(pcm2.Py().native()));
        CHECK(pcm1.Pz().native() == Catch::Approx(pcm2.Pz().native()));
    }
}
