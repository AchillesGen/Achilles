#include "catch2/catch.hpp"
#include "mock_classes.hh"

#include "Achilles/FourVector.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/MomSolver.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Potential.hh"

#include "catch_utils.hh"
#include "spdlog/spdlog.h"

TEST_CASE("CM Momentum Delta solver", "[DeltaFunctions]") {
        achilles::FourVector p1(sqrt(pow(achilles::Constant::mN, 2) + 100*100), 0, 0, 100);
        achilles::FourVector p2(sqrt(pow(achilles::Constant::mN, 2) + 100*100), 0, 0, -100);
        double cosTheta = 0.5;
        double phi = 0;
        auto p3 = achilles::SolveDelta(p1, p2, achilles::Constant::mN, achilles::Constant::mN, cosTheta, phi);
        auto p4 = p1 + p2 - p3;

        double sinTheta = sqrt(0.75);
        CHECK(p1 + p2 == p3 + p4);
        CHECK(p3 == achilles::FourVector(p1.E(), 100*sinTheta, 0, 100*cosTheta));
}

TEST_CASE("Arbitrary Frame Delta solver", "[DeltaFunctions]") {
        std::mt19937 m_rand{std::random_device{}()};
        std::uniform_real_distribution<double> m_dist;

        auto p1 = GENERATE(take(10, randomMomentum(1000, achilles::Constant::mN))); 
        auto p2 = GENERATE(take(10, randomMomentum(1000, achilles::Constant::mN))); 
        double cosTheta = 2*m_dist(m_rand) - 1;
        double phi = 2*M_PI*m_dist(m_rand);

        auto p3 = achilles::SolveDelta(p1, p2, achilles::Constant::mN, achilles::Constant::mN, cosTheta, phi);
        auto p4 = p1 + p2 - p3;
        auto pcm1 = p1+p2;
        auto pcm2 = p3+p4;
        CHECK((pcm1.E() == Approx(pcm2.E()) && pcm1.Px() == Approx(pcm2.Px()) &&
              pcm1.Py() == Approx(pcm2.Py()) && pcm1.Pz() == Approx(pcm2.Pz())));
        CHECK(p3.M() == Approx(achilles::Constant::mN));
        CHECK(p4.M() == Approx(achilles::Constant::mN));
}

TEST_CASE("Potential Delta solver", "[DeltaFunctions]") {
    std::mt19937 m_rand{std::random_device{}()};
    std::uniform_real_distribution<double> m_dist;
    constexpr double rmin = 0.0, rmax = 6.0;
    auto p1 = GENERATE(take(10, randomMomentum(1000, achilles::Constant::mN))); 
    auto p2 = GENERATE(take(10, randomMomentum(225, achilles::Constant::mN))); 
    auto r1 = GENERATE(take(3, random(rmin, rmax)));
    auto r2 = GENERATE(take(3, random(rmin, rmax)));
    auto q_free = p1+p2;
    auto nucleus = std::make_shared<MockNucleus>(); 

    SECTION("WiringaPotential") {
        constexpr double rho0 = 0.16;
        REQUIRE_CALL(*nucleus, Rho(trompeloeil::gt(0)))
            .LR_RETURN((rho0))
            .TIMES(AT_LEAST(4));

        achilles::WiringaPotential potential(rho0);

        auto potential1 = potential(nucleus.get(), p1.P(), r1);
        auto potential2 = potential(nucleus.get(), p2.P(), r2);
        p1.E() = sqrt(p1.P2() + pow(p1.M() + potential1.rscalar, 2)) + potential1.rvector;
        p2.E() = sqrt(p2.P2() + pow(p2.M() + potential2.rscalar, 2)) + potential2.rvector;
        auto q = p1+p2;

        // Rotate so (p1+p2) to be along z-axis
        auto rotation = q.AlignZ();
        q = q.Rotate(rotation);

        // Randomly generate p3Mag and phi
        double phi = 2*M_PI*m_dist(m_rand);
        auto range = achilles::FindMomentumRange(q_free);
        const double dp3 = range.second - range.first;
        achilles::FourVector p3;
        while(true) {
            try {
                double p3Mag = dp3*m_dist(m_rand) + range.first;

                // Solve for p3 and rotate back and calculate p4
                p3 = achilles::SolveDeltaWithPotential(q, nucleus.get(), potential,
                                                       achilles::Constant::mN, achilles::Constant::mN,
                                                       p3Mag, phi, r1, r2);
                break;
            } catch (const std::domain_error &e) {
                continue;
            }
        }
        p3 = p3.RotateBack(rotation);
        auto p4 = p1 + p2 - p3;

        auto pcm1 = p1+p2;
        auto pcm2 = p3+p4;
        CHECK(pcm1.E() == Approx(pcm2.E()));
        CHECK(pcm1.Px() == Approx(pcm2.Px()));
        CHECK(pcm1.Py() == Approx(pcm2.Py()));
        CHECK(pcm1.Pz() == Approx(pcm2.Pz()));
    }

    SECTION("CooperPotential") {
        constexpr size_t AA = 12;
        REQUIRE_CALL(*nucleus, NNucleons())
            .LR_RETURN((AA))
            .TIMES(AT_LEAST(4));

        achilles::CooperPotential potential;

        auto potential1 = potential(nucleus.get(), p1.P(), 1);
        auto potential2 = potential(nucleus.get(), p2.P(), 1);
        p1.E() = sqrt(p1.P2() + pow(p1.M() + potential1.rscalar, 2)) + potential1.rvector;
        p2.E() = sqrt(p2.P2() + pow(p2.M() + potential2.rscalar, 2)) + potential2.rvector;
        auto q = p1+p2;

        // Rotate so (p1+p2) to be along z-axis
        auto rotation = q.AlignZ();
        q = q.Rotate(rotation);

        // Randomly generate p3Mag and phi
        double phi = 2*M_PI*m_dist(m_rand);
        auto range = achilles::FindMomentumRange(q_free);
        const double dp3 = range.second - range.first;
        achilles::FourVector p3;
        while(true) {
            try {
                double p3Mag = dp3*m_dist(m_rand) + range.first;

                // Solve for p3 and rotate back and calculate p4
                p3 = achilles::SolveDeltaWithPotential(q, nucleus.get(), potential, achilles::Constant::mN,
                                                       achilles::Constant::mN, p3Mag, phi, 1, 1);
                break;
            } catch (const std::domain_error &e) {
                continue;
            }
        }
        p3 = p3.RotateBack(rotation);
        auto p4 = p1 + p2 - p3;

        auto pcm1 = p1+p2;
        auto pcm2 = p3+p4;
        CHECK(pcm1.E() == Approx(pcm2.E()));
        CHECK(pcm1.Px() == Approx(pcm2.Px()));
        CHECK(pcm1.Py() == Approx(pcm2.Py()));
        CHECK(pcm1.Pz() == Approx(pcm2.Pz()));
    }
}
