#include "catch2/catch.hpp" 

#include "nuchic/FinalStateMapper.hh"
#include "nuchic/ParticleInfo.hh"
#include "nuchic/FourVector.hh"

TEST_CASE("TwoBodyMapper", "[PhaseSpace]") {
    SECTION("Only works for 2->2") {
        CHECK_THROWS_WITH(nuchic::TwoBodyMapper::Construct({0, 0, 0}),
                          "Incorrect number of masses. Expected 2. Got 3");
        CHECK(nuchic::TwoBodyMapper::Construct({0, 0}) -> NDims() == 2);
    }

    auto mapper = nuchic::TwoBodyMapper::Construct({0, 0});

    SECTION("Forward Map") {
        SECTION("TwoBodyMapper") {
            std::vector<nuchic::FourVector> mom = {{1000, 0, 0, 1000}, {100, 0, 0, -100}, {}, {}};
            std::vector<double> ran = {0.5, 0.5};
            mapper -> GeneratePoint(mom, ran);
            std::vector<double> ran2(2);
            mapper -> GenerateWeight(mom, ran2);
            CHECK(ran[0] == Approx(ran2[0]));
            CHECK(ran[1] == Approx(ran2[1]));
            // TODO: Validate wgt
        }
    }

    SECTION("Reverse Map") {
        SECTION("TwoBodyMapper") {
            const double sqrts = 200;
            std::vector<nuchic::FourVector> mom = {{100, 0, 0, 100},
                                                   {100, 0, 0, -100},
                                                   {sqrts/2, sqrts/2*sin(M_PI/4), 0, sqrts/2*cos(M_PI/4)},
                                                   {sqrts/2, -sqrts/2*sin(M_PI/4), 0, -sqrts/2*cos(M_PI/4)}}; 
            std::vector<double> ran(2);
            mapper -> GenerateWeight(mom, ran);
            std::vector<nuchic::FourVector> mom2(4);
            mom2[0] = mom[0];
            mom2[1] = mom[1];
            mapper -> GeneratePoint(mom2, ran);
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
