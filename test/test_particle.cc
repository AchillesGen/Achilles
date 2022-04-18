#include <sstream>

#include "catch2/catch.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/Particle.hh"
#include "Approx.hh"

constexpr double energy = 1000;
constexpr double timestep = 10000;

TEST_CASE("Formation Zone", "[Particle]") {
    achilles::Particle part{achilles::PID::proton(), {energy, 100, 0, 0}};
    achilles::FourVector mom{energy, 0, 100, 0};

    part.SetFormationZone(part.Momentum(), mom);
    CHECK(part.InFormationZone());
    // Large enough number to ensure it exits the formation zone
    part.UpdateFormationZone(timestep);
    CHECK(!part.InFormationZone());
}

TEST_CASE("Properties", "[Particle]") {
    achilles::Particle part{achilles::PID::proton(), {energy, 100, 0, 0},
                            {0, 1, 0}};

    SECTION("Momentum") {
        CHECK(part.Momentum() == achilles::FourVector{energy, 100, 0, 0});
        part.SetMomentum({energy, 0, 100, 0});
        CHECK(part.Momentum() == achilles::FourVector{energy, 0, 100, 0});
    }

    SECTION("Position") {
        CHECK(part.Position() == achilles::ThreeVector{0, 1, 0});
        part.SetPosition({0, 0, 1});
        CHECK(part.Position() == achilles::ThreeVector{0, 0, 1});
    }

    SECTION("Status") {
        CHECK(part.Status() == achilles::ParticleStatus::background);
        CHECK(part.IsBackground());

        part.Status() = achilles::ParticleStatus::propagating;
        CHECK(part.IsPropagating());

        part.Status() = achilles::ParticleStatus::escaped;
        CHECK(part.IsFinal());
    }

    SECTION("History") {
         
    }
}

TEST_CASE("Propagate", "[Particle]") {
    achilles::Particle part{achilles::PID::proton(), {energy, 100, 0, 0}};
    static constexpr double eps = 1e-8;

    SECTION("Time propagation") {
        part.Propagate(1);
        CHECK_THAT(part.Position(),
                   IsVectorApprox<achilles::ThreeVector>(achilles::ThreeVector{achilles::Constant::HBARC/10, 0, 0}).margin(eps));
        part.BackPropagate(1);
        CHECK_THAT(part.Position(),
                   IsVectorApprox<achilles::ThreeVector>(achilles::ThreeVector{0, 0, 0}).margin(eps));
    }

    SECTION("Space propagation") {
        part.SpacePropagate(1);
        CHECK_THAT(part.Position(),
                   IsVectorApprox<achilles::ThreeVector>(achilles::ThreeVector{1, 0, 0}).margin(eps));
    }
}

TEST_CASE("I/O", "[Particle]") {
    achilles::Particle part{achilles::PID::proton(), {energy, 100, 0, 0}};
    achilles::Particle part2;

    CHECK(part.ToString() == "Particle(2212, FourVector(1000.000000, 100.000000, 0.000000, 0.000000), ThreeVector(0.000000, 0.000000, 0.000000), 0)");

    std::stringstream ss;
    ss << part;
    ss >> part2;

    CHECK(part == part2);
}
