#include "catch2/catch.hpp"

#include "Achilles/FinalStateMapper.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Event.hh"

TEST_CASE("TwoBodyMapper", "[PhaseSpace]") {
    SECTION("Only works for 2->2") {
        CHECK_THROWS_WITH(achilles::TwoBodyMapper::Construct({0, 0, 0}, 1),
                          "Incorrect number of masses. Expected 2. Got 3");
        CHECK(achilles::TwoBodyMapper::Construct({0, 0}, 1)->NDims() == 2);
    }

    auto mapper = achilles::TwoBodyMapper::Construct({0, 0}, 2);

    SECTION("Forward Map") {
        SECTION("TwoBodyMapper") {
            std::vector<achilles::FourVector> mom = {{1000, 0, 0, 1000}, {100, 0, 0, -100}, {}, {}};
			achilles::Event event(mom);
            std::vector<double> ran = {0.5, 0.5};
            mapper->GeneratePoint(event, ran);
            std::vector<double> ran2(2);
            mapper->GenerateWeight(event, ran2);
            CHECK(ran[0] == Approx(ran2[0]));
            CHECK(ran[1] == Approx(ran2[1]));
            // TODO: Validate wgt
        }
    }

    SECTION("Reverse Map") {
        SECTION("TwoBodyMapper") {
            const double sqrts = 200;
            std::vector<achilles::FourVector> mom = {
                {100, 0, 0, 100},
                {100, 0, 0, -100},
                {sqrts / 2, sqrts / 2 * sin(M_PI / 4), 0, sqrts / 2 * cos(M_PI / 4)},
                {sqrts / 2, -sqrts / 2 * sin(M_PI / 4), 0, -sqrts / 2 * cos(M_PI / 4)}};
			achilles::Event event(mom);
            std::vector<double> ran(2);
            mapper->GenerateWeight(event, ran);
			mom=event.Momentum();
            std::vector<achilles::FourVector> mom2(4);
            mom2[0] = mom[0];
            mom2[1] = mom[1];
			achilles::Event event2(mom2);
            mapper->GeneratePoint(event2, ran);
			mom2=event2.Momentum();
            CHECK(mom[2].Px() == Approx(mom2[2].Px()));
            CHECK(mom[2].Py() == Approx(mom2[2].Py()));
            CHECK(mom[2].Pz() == Approx(mom2[2].Pz()));
            CHECK(mom[2].E() == Approx(mom2[2].E()));
            CHECK(mom[3].Px() == Approx(mom2[3].Px()));
            CHECK(mom[3].Py() == Approx(mom2[3].Py()));
            CHECK(mom[3].Pz() == Approx(mom2[3].Pz()));
            CHECK(mom[3].E() == Approx(mom2[3].E()));
        }
    }
}
