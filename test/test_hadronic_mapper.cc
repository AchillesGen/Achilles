#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#include "catch2/catch.hpp"
#pragma GCC diagnostic pop

#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/Event.hh"

TEST_CASE("HadronicMapper", "[PhaseSpace]") {
    achilles::ProcessInfo info({achilles::PID::electron(), {achilles::PID::electron()}});

    SECTION("Forward Map") {
        SECTION("Coherent") {
            info.m_hadronic = {{achilles::PID::carbon()}, {achilles::PID::carbon()}};
            auto mapper = achilles::CoherentMapper::Construct(info, 0);
			std::vector<achilles::FourVector> mom(1);
            achilles::Event event(mom);
            std::vector<double> ran(mapper->NDims());
            mapper->GeneratePoint(event, ran);
            std::vector<double> ran2(mapper->NDims());
            auto wgt = mapper->GenerateWeight(event, ran2);
            CHECK(wgt == 1);
            CHECK(ran == ran2);
        }
        SECTION("QESpectral") {
            info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
            auto mapper = achilles::QESpectralMapper::Construct(info, 1);
			std::vector<achilles::FourVector> mom={{1000, 0, 0, 1000}, {}};
            achilles::Event event(mom);
            std::vector<double> ran = {0.5, 0.5, 0.5, 0.5};
            mapper->SetMasses({0, 0, 0, 0});
            mapper->GeneratePoint(event, ran);
            std::vector<double> ran2(4);
            mapper->GenerateWeight(event, ran2);
            CHECK(ran == ran2);
            // TODO: Validate wgt
        }
    }

    SECTION("Reverse Map") {
        SECTION("Coherent") {
            info.m_hadronic = {{achilles::PID::carbon()}, {achilles::PID::carbon()}};
            auto mapper = achilles::CoherentMapper::Construct(info, 0);
			std::vector<achilles::FourVector> mom={{achilles::ParticleInfo(achilles::PID::carbon()).Mass(), 0, 0, 0}};
            achilles::Event event(mom);
            std::vector<double> ran(mapper->NDims());
            mapper->GenerateWeight(event, ran);
			std::vector<achilles::FourVector> mom2(1);
            achilles::Event event2(mom2);
            mapper->GeneratePoint(event2, ran);
            CHECK(event == event2);
        }
        SECTION("QESpectral") {
            info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
            auto mapper = achilles::QESpectralMapper::Construct(info, 1);
			std::vector<achilles::FourVector> mom={
                {1000, 0, 0, 1000},
                {achilles::Constant::mN - 20, 400 * sin(M_PI / 4) * cos(M_PI / 6),
                 400 * sin(M_PI / 4) * sin(M_PI / 6), 400 * cos(M_PI / 4)},
            };
            achilles::Event event(mom);
            std::vector<double> ran(4);
            mapper->SetMasses({0, 0, 0, 0});
            mapper->GenerateWeight(event, ran);
			std::vector<achilles::FourVector> mom2(2);
            achilles::Event event2(mom2);
            event2.Momentum()[0] = event.Momentum()[0];
            mapper->GeneratePoint(event2, ran);
            CHECK(event.Momentum()[1].Px() == Approx(event2.Momentum()[1].Px()));
            CHECK(event.Momentum()[1].Py() == Approx(event2.Momentum()[1].Py()));
            CHECK(event.Momentum()[1].Pz() == Approx(event2.Momentum()[1].Pz()));
            CHECK(event.Momentum()[1].E() == Approx(event2.Momentum()[1].E()));
        }
    }
}
