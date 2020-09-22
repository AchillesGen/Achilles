#include "catch2/catch.hpp"

#include "nuchic/FourVector.hh"
#include "nuchic/Rambo.hh"
#include "nuchic/Random.hh"

TEST_CASE("Phase Space Momentums", "[Rambo]") {
    SECTION("Massless Two Body") {
        nuchic::Rambo rambo(2);
        randutils::mt19937_rng rng{};

        std::vector<double> rans(8);
        rng.generate(rans, 0.0, 1.0);

        std::vector<nuchic::FourVector> momentums = {
            {0, 0, 8000, 8000}, {0, 0, -8000, 8000},
            {0, 0, 0, 0}, {0, 0, 0, 0}
        };

        std::vector<double> masses = {0, 0, 0, 0};

        rambo.GeneratePoint(momentums, masses, rans);

        CHECK(momentums[2].Px() == Approx(-momentums[3].Px()).margin(1e-8));
        CHECK(momentums[2].Py() == Approx(-momentums[3].Py()).margin(1e-8));
        CHECK(momentums[2].Pz() == Approx(-momentums[3].Pz()).margin(1e-8));
        CHECK(momentums[2].E() == Approx(momentums[3].E()).margin(1e-8));

        auto momIn = momentums[0] + momentums[1];
        auto momOut = momentums[2] + momentums[3];
        CHECK(momIn.Px() == Approx(momOut.Px()).margin(1e-8));
        CHECK(momIn.Py() == Approx(momOut.Py()).margin(1e-8));
        CHECK(momIn.Pz() == Approx(momOut.Pz()).margin(1e-8));
        CHECK(momIn.E() == Approx(momOut.E()).margin(1e-8));
    }

    SECTION("Massive Two Body") {
        nuchic::Rambo rambo(2);
        randutils::mt19937_rng rng{};

        std::vector<double> rans(8);
        rng.generate(rans, 0.0, 1.0);

        std::vector<nuchic::FourVector> momentums = {
            {0, 0, 8000, 8000}, {0, 0, -8000, 8000},
            {0, 0, 0, 0}, {0, 0, 0, 0}
        };

        std::vector<double> masses(4);
        rng.generate(masses, 0.0, 2000.0);

        rambo.GeneratePoint(momentums, masses, rans);

        CHECK(momentums[2].Px() == Approx(-momentums[3].Px()).margin(1e-8));
        CHECK(momentums[2].Py() == Approx(-momentums[3].Py()).margin(1e-8));
        CHECK(momentums[2].Pz() == Approx(-momentums[3].Pz()).margin(1e-8));
        CHECK(momentums[2].E()*momentums[2].E()-masses[2]*masses[2]
              == Approx(momentums[3].E()*momentums[3].E()-masses[3]*masses[3]).margin(1e-8));

        auto momIn = momentums[0] + momentums[1];
        auto momOut = momentums[2] + momentums[3];
        CHECK(momIn.Px() == Approx(momOut.Px()).margin(1e-8));
        CHECK(momIn.Py() == Approx(momOut.Py()).margin(1e-8));
        CHECK(momIn.Pz() == Approx(momOut.Pz()).margin(1e-8));
        CHECK(momIn.E() == Approx(momOut.E()).margin(1e-8));
    }

    SECTION("Massless Three Body") {
        nuchic::Rambo rambo(3);
        randutils::mt19937_rng rng{};

        std::vector<double> rans(12);
        rng.generate(rans, 0.0, 1.0);

        std::vector<nuchic::FourVector> momentums = {
            {0, 0, 8000, 8000}, {0, 0, -8000, 8000},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}
        };

        std::vector<double> masses = {0, 0, 0, 0, 0};

        rambo.GeneratePoint(momentums, masses, rans);

        auto momOut1 = momentums[2]+momentums[3];
        CHECK(momOut1.Px() == Approx(-momentums[4].Px()).margin(1e-8));
        CHECK(momOut1.Py() == Approx(-momentums[4].Py()).margin(1e-8));
        CHECK(momOut1.Pz() == Approx(-momentums[4].Pz()).margin(1e-8));
        CHECK(momOut1.E() + momentums[4].E() == Approx(16000).margin(1e-8));

        auto momIn = momentums[0] + momentums[1];
        auto momOut = momentums[2] + momentums[3] + momentums[4];
        CHECK(momIn.Px() == Approx(momOut.Px()).margin(1e-8));
        CHECK(momIn.Py() == Approx(momOut.Py()).margin(1e-8));
        CHECK(momIn.Pz() == Approx(momOut.Pz()).margin(1e-8));
        CHECK(momIn.E() == Approx(momOut.E()).margin(1e-8));
    }

    SECTION("Massive Three Body") {
        nuchic::Rambo rambo(3);
        randutils::mt19937_rng rng{};

        std::vector<double> rans(12);
        rng.generate(rans, 0.0, 1.0);

        std::vector<nuchic::FourVector> momentums = {
            {0, 0, 8000, 8000}, {0, 0, -8000, 8000},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}
        };

        std::vector<double> masses(momentums.size());
        rng.generate(masses, 0.0, 2000.0);

        rambo.GeneratePoint(momentums, masses, rans);

        auto momOut1 = momentums[2]+momentums[3];
        CHECK(momOut1.Px() == Approx(-momentums[4].Px()).margin(1e-8));
        CHECK(momOut1.Py() == Approx(-momentums[4].Py()).margin(1e-8));
        CHECK(momOut1.Pz() == Approx(-momentums[4].Pz()).margin(1e-8));
        CHECK(momOut1.E() + momentums[4].E() == Approx(16000).margin(1e-8));

        auto momIn = momentums[0] + momentums[1];
        auto momOut = momentums[2] + momentums[3] + momentums[4];
        CHECK(momIn.Px() == Approx(momOut.Px()).margin(1e-8));
        CHECK(momIn.Py() == Approx(momOut.Py()).margin(1e-8));
        CHECK(momIn.Pz() == Approx(momOut.Pz()).margin(1e-8));
        CHECK(momIn.E() == Approx(momOut.E()).margin(1e-8));
    }
}

TEST_CASE("Phase Space Weights", "[Rambo]") {
    SECTION("Massless Two Body") {
        nuchic::Rambo rambo(2);
        randutils::mt19937_rng rng{};

        std::vector<double> rans(8);
        rng.generate(rans, 0.0, 1.0);

        std::vector<nuchic::FourVector> momentums = {
            {0, 0, 8000, 8000}, {0, 0, -8000, 8000},
            {0, 0, 0, 0}, {0, 0, 0, 0}
        };

        std::vector<double> masses = {0, 0, 0, 0};

        rambo.GeneratePoint(momentums, masses, rans);
        rambo.GenerateWeight(momentums, masses);

        double weight = M_PI/2/pow(2*M_PI, 2);
        CHECK(rambo.Weight() == Approx(weight).margin(1e-8));
    }

    SECTION("Massive Two Body") {
        nuchic::Rambo rambo(2);
        randutils::mt19937_rng rng{};

        std::vector<double> rans(8);
        rng.generate(rans, 0.0, 1.0);

        std::vector<nuchic::FourVector> momentums = {
            {0, 0, 8000, 8000}, {0, 0, -8000, 8000},
            {0, 0, 0, 0}, {0, 0, 0, 0}
        };

        std::vector<double> masses = {1000, 1000, 5000, 100};

        rambo.GeneratePoint(momentums, masses, rans);
        rambo.GenerateWeight(momentums, masses);

        double weight = M_PI/2/pow(2*M_PI, 2);
        CHECK(rambo.Weight() == Approx(weight).margin(1e-8));
    }

    SECTION("Massless Three Body") {
        nuchic::Rambo rambo(3);
        randutils::mt19937_rng rng{};

        std::vector<double> rans(12);
        rng.generate(rans, 0.0, 1.0);

        std::vector<nuchic::FourVector> momentums = {
            {0, 0, 8000, 8000}, {0, 0, -8000, 8000},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}
        };

        std::vector<double> masses = {0, 0, 0, 0, 0};

        rambo.GeneratePoint(momentums, masses, rans);
        rambo.GenerateWeight(momentums, masses);

        double weight = exp(2*log(16000)+2*log(M_PI/2)-log(2))/pow(2*M_PI, 5);
        CHECK(rambo.Weight() == Approx(weight).margin(1e-8));
    }
}
