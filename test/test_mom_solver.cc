#include "catch2/catch.hpp"

#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"
#include "nuchic/MomSolver.hh"
#include "nuchic/Constants.hh"

#include "catch_utils.hh"
#include "spdlog/spdlog.h"

TEST_CASE("Delta function solver", "[Cascade]") {
    SECTION("Works in the CM") {
        nuchic::FourVector p1(0, 0, 100, sqrt(pow(nuchic::Constant::mN, 2) + 100*100));
        nuchic::FourVector p2(0, 0, -100, sqrt(pow(nuchic::Constant::mN, 2) + 100*100));
        double cosTheta = 0.5;
        double phi = 0;
        auto p3 = nuchic::SolveDelta(p1, p2, nuchic::Constant::mN, nuchic::Constant::mN, cosTheta, phi);
        auto p4 = p1 + p2 - p3;

        double sinTheta = sqrt(0.75);
        CHECK(p1 + p2 == p3 + p4);
        CHECK(p3 == nuchic::FourVector(100*sinTheta, 0, 100*cosTheta, p1.E()));
    }

    SECTION("Works for random momentum") {
        std::mt19937 m_rand{std::random_device{}()};
        std::uniform_real_distribution<double> m_dist;

        auto p1 = GENERATE(take(1, randomMomentum(1000, nuchic::Constant::mN))); 
        auto p2 = GENERATE(take(1, randomMomentum(1000, nuchic::Constant::mN))); 
        double cosTheta = 2*m_dist(m_rand) - 1;
        double phi = 2*M_PI*m_dist(m_rand);

        auto p3 = nuchic::SolveDelta(p1, p2, nuchic::Constant::mN, nuchic::Constant::mN, cosTheta, phi);
        auto p4 = p1 + p2 - p3;
        auto pcm1 = p1+p2;
        auto pcm2 = p3+p4;
        CHECK((pcm1.E() == Approx(pcm2.E()) && pcm1.Px() == Approx(pcm2.Px()) &&
              pcm1.Py() == Approx(pcm2.Py()) && pcm1.Pz() == Approx(pcm2.Pz())));
        CHECK(p3.M() == Approx(nuchic::Constant::mN));
        CHECK(p4.M() == Approx(nuchic::Constant::mN));
    }

    SECTION("Works with random momentum and Potentials") {
        std::mt19937 m_rand{std::random_device{}()};
        std::uniform_real_distribution<double> m_dist;

        auto p1 = GENERATE(take(100, randomMomentum(1000, nuchic::Constant::mN))); 
        auto p2 = GENERATE(take(100, randomMomentum(225, nuchic::Constant::mN))); 
        auto q_free = p1+p2;
        p1.E() += nuchic::Potential(p1.P(), 0.16); 
        p2.E() += nuchic::Potential(p2.P(), 0.16); 
        auto q = p1+p2;

        // Rotate so (p1+p2) to be along z-axis
        auto rotation = q.AlignZ();
        q = q.Rotate(rotation);

        // Randomly generate p3Mag and phi
        double phi = 2*M_PI*m_dist(m_rand);
        std::pair<double, double> range;
        const double detM = sqrt((q_free.M2()-4*pow(nuchic::Constant::mN, 2))*(q_free.M2() + q_free.P2())/q_free.M2());
        range.first = q_free.M2() > 2*nuchic::Constant::mN*q_free.E() ? (detM-q_free.P())/2 : (q_free.P()-detM)/2;
        range.second = (q_free.P()+detM)/2;
        constexpr double rangeExtend = 1.05;
        range.first /= rangeExtend;
        range.second *= rangeExtend;
        const double dp3 = range.second - range.first;
        nuchic::FourVector p3;
        size_t iteration = 0;
        static size_t iteration_sum = 0;
        while(true) {
            spdlog::info("Iteration: {}", iteration++);
            try {
                double p3Mag = dp3*m_dist(m_rand) + range.first;
                spdlog::info("p3Mag = {}", p3Mag);
                if(iteration != 1) {
                    spdlog::info("q = {}", q);
                    spdlog::info("q_free = {}", q_free);
                    spdlog::info("p3Range = [{}, {}]", range.first, range.second);
                    spdlog::info("q = {}", q_free.P());
                    spdlog::info("M = {}", q_free.M());
                    spdlog::info("U = {}", nuchic::Potential(p1.P(), 0.16) + nuchic::Potential(p2.P(), 0.16));
                }
                // spdlog::info("p3Mag = {}", p3Mag);
                // const double cosGuess = -(q.M2()-2*q.E()*sqrt(p3Mag*p3Mag+pow(nuchic::Constant::mN, 2)))/(2*q.P()*p3Mag);
                // spdlog::info("cosTheta_guess = {}", cosGuess);

                // Solve for p3 and rotate back and calculate p4
                p3 = nuchic::SolveDeltaWithPotential(q, nuchic::Constant::mN, nuchic::Constant::mN,
                                                     p3Mag, phi, 0.16, 0.16);
                break;
            } catch (const std::domain_error &e) {
                continue;
            }
        }
        iteration_sum += iteration;
        p3 = p3.RotateBack(rotation);
        auto p4 = p1 + p2 - p3;

        spdlog::debug("p1 = {}", p1);
        spdlog::debug("p2 = {}", p2);
        spdlog::debug("p3 = {}", p3);
        spdlog::debug("p4 = {}", p4);

        auto pcm1 = p1+p2;
        auto pcm2 = p3+p4;
        CHECK(pcm1.E() == Approx(pcm2.E()));
        CHECK(pcm1.Px() == Approx(pcm2.Px()));
        CHECK(pcm1.Py() == Approx(pcm2.Py()));
        CHECK(pcm1.Pz() == Approx(pcm2.Pz()));

        spdlog::info("Total iterations: {}", iteration_sum);
    }
}
