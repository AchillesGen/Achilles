#include <sstream>

#include "catch2/catch.hpp"

#include "nuchic/Particle.hh"
#include "Approx.hh"

TEST_CASE("Formation Zone", "[Particle]") {
    nuchic::Particle part{2212, {100, 0, 0, 1000}};
    nuchic::FourVector mom{0, 100, 0, 1000};

    part.SetFormationZone(part.Momentum(), mom);
    CHECK(part.InFormationZone());
    // Large enough number to ensure it exits the formation zone
    part.UpdateFormationZone(10000);
    CHECK(!part.InFormationZone());
}

TEST_CASE("Properties", "[Particle]") {
    nuchic::Particle part{nuchic::PID::proton(), {100, 0, 0, 1000},
                          {0, 1, 0}};

    SECTION("Momentum") {
        CHECK(part.Momentum() == nuchic::FourVector{100, 0, 0, 1000});
        part.SetMomentum({0, 100, 0, 1000});
        CHECK(part.Momentum() == nuchic::FourVector{0, 100, 0, 1000});
    }

    SECTION("Position") {
        CHECK(part.Position() == nuchic::ThreeVector{0, 1, 0});
        part.SetPosition({0, 0, 1});
        CHECK(part.Position() == nuchic::ThreeVector{0, 0, 1});
    }

    SECTION("Status") {
        CHECK(part.Status() == nuchic::ParticleStatus::background);
        CHECK(part.IsBackground());

        part.Status() = nuchic::ParticleStatus::propagating;
        CHECK(part.IsPropagating());

        part.Status() = nuchic::ParticleStatus::escaped;
        CHECK(part.IsFinal());
    }

    SECTION("History") {
         
    }
}

TEST_CASE("Propagate", "[Particle]") {
    nuchic::Particle part{nuchic::PID::proton(), {100, 0, 0, 1000}};
    static constexpr double eps = 1e-8;

    SECTION("Time propagation") {
        part.Propagate(1);
        CHECK_THAT(part.Position(),
                   IsVectorApprox<nuchic::ThreeVector>(nuchic::ThreeVector{20, 0, 0}).margin(eps));
        part.BackPropagate(1);
        CHECK_THAT(part.Position(),
                   IsVectorApprox<nuchic::ThreeVector>(nuchic::ThreeVector{0, 0, 0}).margin(eps));
    }

    SECTION("Space propagation") {
        part.SpacePropagate(1);
        CHECK_THAT(part.Position(),
                   IsVectorApprox<nuchic::ThreeVector>(nuchic::ThreeVector{1, 0, 0}).margin(eps));
    }
}

TEST_CASE("I/O", "[Particle]") {
    nuchic::Particle part{nuchic::PID::proton(), {100, 0, 0, 1000}};
    nuchic::Particle part2;

    CHECK(part.ToString() == "Particle(2212, FourVector(100.000000, 0.000000, 0.000000, 1000.000000), ThreeVector(0.000000, 0.000000, 0.000000), 0)");

    std::stringstream ss;
    ss << part;
    ss >> part2;

    CHECK(part == part2);
}
