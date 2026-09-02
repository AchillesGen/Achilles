// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/ProcessInfo.hh"

namespace units = achilles::units;
using namespace achilles::units::literals;

TEST_CASE("HadronicMapper", "[PhaseSpace]") {
    achilles::ProcessInfo info({achilles::PID::electron(), {achilles::PID::electron()}});

    SECTION("Forward Map") {
        SECTION("Coherent") {
            info.m_hadronic = {{achilles::PID::carbon()}, {achilles::PID::carbon()}};
            auto mapper = achilles::CoherentMapper::Construct(info, 0);
            std::vector<achilles::FourVector> mom(1);
            std::vector<double> ran(mapper->NDims());
            mapper->GeneratePoint(mom, ran);
            std::vector<double> ran2(mapper->NDims());
            auto wgt = mapper->GenerateWeight(mom, ran2);
            CHECK(wgt == 1);
            CHECK(ran == ran2);
        }
        SECTION("QESpectral") {
            info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
            auto mapper = achilles::QESpectralMapper::Construct(info, 1);
            std::vector<achilles::FourVector> mom = {{1000_MeV, 0_MeV, 0_MeV, 1000_MeV}, {}};
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
            info.m_hadronic = {{achilles::PID::carbon()}, {achilles::PID::carbon()}};
            auto mapper = achilles::CoherentMapper::Construct(info, 0);
            std::vector<achilles::FourVector> mom = {
                {units::Energy{achilles::ParticleInfo(achilles::PID::carbon()).Mass().native()},
                 0_MeV, 0_MeV, 0_MeV}};
            std::vector<double> ran(mapper->NDims());
            mapper->GenerateWeight(mom, ran);
            std::vector<achilles::FourVector> mom2(1);
            mapper->GeneratePoint(mom2, ran);
            CHECK(mom == mom2);
        }
        SECTION("QESpectral") {
            info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
            auto mapper = achilles::QESpectralMapper::Construct(info, 1);
            std::vector<achilles::FourVector> mom = {
                {1000_MeV, 0_MeV, 0_MeV, 1000_MeV},
                {achilles::Constant::mN - 20_MeV, 400 * sin(M_PI / 4) * cos(M_PI / 6) * units::MeV,
                 400 * sin(M_PI / 4) * sin(M_PI / 6) * units::MeV,
                 400 * cos(M_PI / 4) * units::MeV},
            };
            std::vector<double> ran(4);
            mapper->SetMasses({0, 0, 0, 0});
            mapper->GenerateWeight(mom, ran);
            std::vector<achilles::FourVector> mom2(2);
            mom2[0] = mom[0];
            mapper->GeneratePoint(mom2, ran);
            CHECK(mom[1].Px().native() == Catch::Approx(mom2[1].Px().native()));
            CHECK(mom[1].Py().native() == Catch::Approx(mom2[1].Py().native()));
            CHECK(mom[1].Pz().native() == Catch::Approx(mom2[1].Pz().native()));
            CHECK(mom[1].E().native() == Catch::Approx(mom2[1].E().native()));
        }
    }
}
