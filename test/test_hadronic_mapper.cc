#include "catch2/catch.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/ParticleInfo.hh"

TEST_CASE("HadronicMapper", "[PhaseSpace]") {
    const double mCarbon = achilles::ParticleInfo(achilles::PID::carbon()).Mass();
    SECTION("Forward Map") {
        SECTION("Coherent") {
            auto mapper = achilles::CoherentMapper::Construct(0);
            mapper->SetMasses({pow(mCarbon, 2)});
            std::vector<achilles::FourVector> mom(1);
            std::vector<double> ran(mapper->NDims());
            mapper->GeneratePoint(mom, ran);
            std::vector<double> ran2(mapper->NDims());
            auto wgt = mapper->GenerateWeight(mom, ran2);
            CHECK(wgt == 1);
            CHECK(ran == ran2);
        }
        SECTION("QESpectral") {
            auto mapper = achilles::QESpectralMapper::Construct(1);
            std::vector<achilles::FourVector> mom = {{1000, 0, 0, 1000}, {}};
            std::vector<double> ran = {0.5, 0.5, 0.5, 0.5};
            mapper->SetMasses({0, 0, 0, 0});
            mapper->GeneratePoint(mom, ran);
            std::vector<double> ran2(4);
            mapper->GenerateWeight(mom, ran2);
            CHECK(ran == ran2);
            // TODO: Validate wgt
        }
    }

    SECTION("Reverse Map") {
        SECTION("Coherent") {
            auto mapper = achilles::CoherentMapper::Construct(0);
            mapper->SetMasses({pow(mCarbon, 2)});
            std::vector<achilles::FourVector> mom = {{mCarbon, 0, 0, 0}};
            std::vector<double> ran(mapper->NDims());
            mapper->GenerateWeight(mom, ran);
            std::vector<achilles::FourVector> mom2(1);
            mapper->GeneratePoint(mom2, ran);
            CHECK(mom == mom2);
        }
        SECTION("QESpectral") {
            auto mapper = achilles::QESpectralMapper::Construct(1);
            std::vector<achilles::FourVector> mom = {
                {1000, 0, 0, 1000},
                {achilles::Constant::mN - 20, 400 * sin(M_PI / 4) * cos(M_PI / 6),
                 400 * sin(M_PI / 4) * sin(M_PI / 6), 400 * cos(M_PI / 4)},
            };
            std::vector<double> ran(4);
            mapper->SetMasses({0, 0, 0, 0});
            mapper->GenerateWeight(mom, ran);
            std::vector<achilles::FourVector> mom2(2);
            mom2[0] = mom[0];
            mapper->GeneratePoint(mom2, ran);
            CHECK(mom[1].Px() == Approx(mom2[1].Px()));
            CHECK(mom[1].Py() == Approx(mom2[1].Py()));
            CHECK(mom[1].Pz() == Approx(mom2[1].Pz()));
            CHECK(mom[1].E() == Approx(mom2[1].E()));
        }
    }
}
