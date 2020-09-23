#include "catch2/catch.hpp"

#include <sstream>

#include "nuchic/Particle.hh"

TEST_CASE("Particle Initialization", "[Particle]") {
    SECTION("Constructors") {
        nuchic::Particle p;

        CHECK(p.ID() == nuchic::PID::undefined());
        CHECK(p.Momentum() == nuchic::FourVector());
        CHECK(p.Status() == nuchic::ParticleStatus::background);
    }

    SECTION("Assignment") {
        nuchic::Particle p;

        p = nuchic::Particle(nuchic::PID::neutron());
        CHECK(p.ID() == nuchic::PID::neutron());

        nuchic::Particle p2(p);
        CHECK(p2 == p);
    }
}

TEST_CASE("Particle Setters and Getters", "[Particle]") {
    SECTION("Trivial Setters and Getters") {
        nuchic::Particle p;
        p.SetPID(nuchic::PID::proton());
        p.SetPosition({10, 10, 10});
        p.SetMomentum({10, 0, 0, 1100});
        p.SetStatus(nuchic::ParticleStatus::captured);
        p.SetMothers({0, 1, 2});
        p.SetDaughters({3, 4});

        CHECK(p.ID() == nuchic::PID::proton());
        CHECK(p.Position() == nuchic::ThreeVector(10, 10, 10));
        CHECK(p.Momentum() == nuchic::FourVector(10, 0, 0, 1100));
        CHECK(p.Status() == nuchic::ParticleStatus::captured);
        CHECK(p.Mothers() == std::vector<int>({0, 1, 2}));
        CHECK(p.Daughters() == std::vector<int>({3, 4}));

        p.AddMother(5);
        p.AddDaughter(6);
        CHECK(p.Mothers() == std::vector<int>({0, 1, 2, 5}));
        CHECK(p.Daughters() == std::vector<int>({3, 4, 6}));

        CHECK(!p.IsBackground());
        CHECK(!p.IsPropagating());
        CHECK(!p.IsFinal());
        CHECK(p.GetDistanceTraveled() == 0);

        CHECK(p.Info().ID() == p.ID());
    }

    SECTION("Boost") {
        nuchic::Particle p(2212, {100, 0, 0, 1100});

        p.SetMomentum(p.Momentum().Boost(-p.Beta()));
        CHECK(p.Momentum() == nuchic::FourVector(0, 0, 0, p.Mass()));
        CHECK(p.Px() == 0);
        CHECK(p.Py() == 0);
        CHECK(p.Pz() == 0);
        CHECK(p.E() == p.Mass());
    }
}

TEST_CASE("Formation Zone", "[Particle]") {
    nuchic::Particle p0(nuchic::PID::proton(), {100, 0, 0, 1100});
    nuchic::Particle p1(nuchic::PID::proton(), {0, 100, 0, 1100});
    CHECK(p1 != p0);

    p1.SetFormationZone(p0.Momentum(), p1.Momentum());
    CHECK(p1.InFormationZone());
    CHECK(p1.FormationZone() > 0);
    p1.UpdateFormationZone(p1.FormationZone());
    CHECK(!p1.InFormationZone());
}

TEST_CASE("Propagate", "[Particle]") {
    SECTION("Time Propagate") {
        nuchic::Particle p(nuchic::PID::proton(), {100, 0, 0, 1100});

        auto position0 = p.Position();
        p.Propagate(0.2);
        auto radius = p.Radius();
        CHECK(radius > 0);
        p.Propagate(0.2);
        CHECK(radius < p.Radius());
        p.BackPropagate(0.4);
        CHECK(p.Position() == position0);
    }

    SECTION("Space Propagate") {
        nuchic::Particle p(nuchic::PID::proton(), {100, 0, 0, 1100});

        auto position0 = p.Position();
        p.SpacePropagate(0.1);
        p.SpacePropagate(0.1);
        p.SpacePropagate(-0.2);
        CHECK(p.Position() == position0);
    }
}

TEST_CASE("Particle I/O", "[Particle]") {
    nuchic::Particle p(nuchic::PID::proton(), {100, 0, 0, 1100});

    std::string test = "Particle(" + std::to_string(int(nuchic::PID::proton())) + ", "
        + nuchic::FourVector(100, 0, 0, 1100).ToString() + ", "
        + nuchic::ThreeVector(0, 0, 0).ToString() + ", " 
        + std::to_string(static_cast<int>(nuchic::ParticleStatus::background)) + ")";

    CHECK(p.ToString() == test);

    std::stringstream ss;
    ss << p;
    nuchic::Particle p2;
    ss >> p2;
    CHECK(p == p2);

}
