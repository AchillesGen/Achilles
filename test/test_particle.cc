// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include <sstream>

#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/Particle.hh"
#include "Approx.hh"

using achilles::FourVector;
using achilles::Particle;
using achilles::ThreePosition;
namespace units = achilles::units;
using namespace achilles::units::literals;

constexpr units::Energy energy = 1000_MeV;
constexpr units::Time timestep = 10000.0 * units::fm;

TEST_CASE("Formation Zone", "[Particle]") {
    Particle part{achilles::PID::proton(), {units::Energy{energy}, 100_MeV, 0_MeV, 0_MeV}};
    FourVector mom{energy, 0_MeV, 100_MeV, 0_MeV};

    part.SetFormationZone(part.Momentum(), mom);
    CHECK(part.InFormationZone());
    // Large enough number to ensure it exits the formation zone
    part.UpdateFormationZone(timestep);
    CHECK(!part.InFormationZone());
}

TEST_CASE("Properties", "[Particle]") {
    Particle part{achilles::PID::proton(),
                  {units::Energy{energy}, 100_MeV, 0_MeV, 0_MeV},
                  {0.0_fm, 1.0_fm, 0.0_fm}};

    SECTION("Momentum") {
        CHECK(part.Momentum() == FourVector{energy, 100_MeV, 0_MeV, 0_MeV});
        part.SetMomentum({energy, 0_MeV, 100_MeV, 0_MeV});
        CHECK(part.Momentum() == FourVector{energy, 0_MeV, 100_MeV, 0_MeV});
    }

    SECTION("Position") {
        CHECK(part.Position() == ThreePosition{0.0_fm, 1.0_fm, 0.0_fm});
        part.SetPosition({0.0_fm, 0.0_fm, 1.0_fm});
        CHECK(part.Position() == ThreePosition{0.0_fm, 0.0_fm, 1.0_fm});
    }

    SECTION("Status") {
        CHECK(part.Status() == achilles::ParticleStatus::background);
        CHECK(part.IsBackground());

        part.Status() = achilles::ParticleStatus::propagating;
        CHECK(part.IsPropagating());

        part.Status() = achilles::ParticleStatus::final_state;
        CHECK(part.IsFinal());
    }

    SECTION("History") {}
}

TEST_CASE("Propagate", "[Particle]") {
    Particle part{achilles::PID::proton(), {units::Energy{energy}, 100_MeV, 0_MeV, 0_MeV}};
    static constexpr double eps = 1e-8;

    SECTION("Time propagation") {
        part.Propagate(1.0 * units::fm);
        const double beta = (part.Momentum().P() / part.E()).native();
        CHECK_THAT(part.Position(),
                   IsVectorApprox<ThreePosition>(ThreePosition{beta * units::fm, 0.0_fm, 0.0_fm})
                       .margin(eps));
        part.BackPropagate(1.0 * units::fm);
        CHECK_THAT(
            part.Position(),
            IsVectorApprox<ThreePosition>(ThreePosition{0.0_fm, 0.0_fm, 0.0_fm}).margin(eps));
    }

    SECTION("Space propagation") {
        part.SpacePropagate(1.0 * units::fm);
        CHECK_THAT(
            part.Position(),
            IsVectorApprox<ThreePosition>(ThreePosition{1.0_fm, 0.0_fm, 0.0_fm}).margin(eps));
    }
}

TEST_CASE("I/O", "[Particle]") {
    Particle part{achilles::PID::proton(), {units::Energy{energy}, 100_MeV, 0_MeV, 0_MeV}};
    Particle part2;

    CHECK(part.ToString() == "Particle(2212, FourVector(1000.000000, 100.000000, 0.000000, "
                             "0.000000), ThreeVector(0.000000, 0.000000, 0.000000), 25)");

    std::stringstream ss;
    ss << part;
    ss >> part2;

    CHECK(part == part2);
}

TEST_CASE("Closest Approach", "[Particle]") {
    // A particle moving at beta = 1/hbarc in z, one fm-scale offset away.
    Particle part1(achilles::PID::proton(),
                   {1_MeV, 0_MeV, 0_MeV, units::Energy{1.0 / achilles::Constant::HBARC}});
    Particle part2(achilles::PID::proton(), {}, {3.0_fm, 2.0_fm, 1.0_fm});

    auto time = achilles::ClosestApproach(part1, part2);
    CHECK(time.in(units::fm) == Catch::Approx(achilles::Constant::HBARC));

    part1.Propagate(time);
    CHECK(part1.Position() == ThreePosition{0.0_fm, 0.0_fm, 1.0_fm});
}
