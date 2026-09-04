// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"

#include "Achilles/ParticleInfo.hh"

TEST_CASE("ParticleInfo", "[ParticleInfo]") {
    SECTION("Must be a valid particle") {
        CHECK_THROWS_WITH(achilles::ParticleInfo(23413), "Invalid PID: id=23413");

        CHECK_THROWS_WITH(achilles::ParticleInfo(achilles::PID(23413)), "Invalid PID: id=23413");

        for(const auto &entry : achilles::ParticleInfo::Database()) {
            CHECK_NOTHROW(achilles::ParticleInfo(entry.first));
        }
    }

    SECTION("Add Entry to database") {
        const auto size = achilles::ParticleInfo::Database().size();
        auto entry = std::make_shared<achilles::ParticleInfoEntry>(
            achilles::PID(123456789), 0, 0, 0, 0, 0, 0, 0, 0, true, false, "test", "anti-test");
        achilles::ParticleInfo info(entry);

        CHECK(achilles::ParticleInfo::Database().size() == size + 1);
    }

    SECTION("Negative PIDs are interpreted as anti-particles") {
        achilles::ParticleInfo info1(achilles::PID::proton()), info2(-2212);

        CHECK(info1.Anti() == info2);
        CHECK(info1.Anti().Anti() == info1);
    }

    SECTION("Valid properties are returned") {
        // dummy particle
        auto entry = std::make_shared<achilles::ParticleInfoEntry>(
            achilles::PID(123456789), 10, 1, 1, 1, 1, 1, 1, 0, true, false, "test", "anti-test");
        achilles::ParticleInfo info(entry);
        achilles::ParticleInfo ainfo(entry, true);

        CHECK(info.Name() == "test");
        CHECK(ainfo.Name() == "anti-test");

        CHECK(info.ID() == achilles::PID(123456789));
        CHECK(ainfo.ID() == achilles::PID(123456789));

        CHECK(info.IntID() == 123456789);
        CHECK(ainfo.IntID() == -123456789);

        CHECK(info.IsBaryon() == true);
        CHECK(info.IsHadron() == false);
        CHECK(info.IsBHadron() == false);
        CHECK(info.IsCHadron() == false);
        CHECK(info.IsAnti() == false);
        CHECK(ainfo.IsAnti() == true);
        CHECK(info.IsFermion() == true);
        CHECK(info.IsBoson() == false);
        CHECK(info.IsScalar() == false);
        CHECK(info.IsVector() == false);
        CHECK(info.IsTensor() == false);
        CHECK(info.IsPhoton() == false);
        CHECK(info.IsLepton() == false);
        CHECK(info.IsQuark() == false);
        CHECK(info.IsGluon() == false);

        CHECK(info.IntCharge() == 1);
        CHECK(ainfo.IntCharge() == -1);
        CHECK(info.Charge() == 1. / 3);
        CHECK(ainfo.Charge() == -1. / 3);
        CHECK(info.IntSpin() == 1);
        CHECK(info.Spin() == 1. / 2);

        CHECK(info.SelfAnti() == false);
        CHECK(info.Majorana() == false);
        CHECK(info.Stable() == 1);
        CHECK(info.IsStable() == true);
        CHECK(info.IsMassive() == true);
        CHECK(info.Mass().native() == 10);
        CHECK(info.Width().native() == 1);
    }

    SECTION("Equality of two particles") {
        achilles::ParticleInfo info1(achilles::PID::proton()), info2(achilles::PID::proton()),
            info3(achilles::PID::neutron()), info4(-2212);

        // Only same PIDs are equal
        CHECK(info1 == info2);
        CHECK(info1 != info3);
        // Anti-particles are not equal to particles
        CHECK(info1 != info4);
    }
}

TEST_CASE("Mass constants come from the particle database", "[ParticleInfo]") {
    using achilles::ParticleInfo;
    using achilles::PID;
    namespace Constant = achilles::Constant;

    // Bit-exact: these were hard-coded constants before they moved into the
    // particle file, and the Fortran interference benchmark holds to 1e-6.
    CHECK(Constant::mp().native() == 938.27208816);
    CHECK(Constant::mn().native() == 939.56542054);

    CHECK(Constant::mp() == ParticleInfo(PID::proton()).Mass());
    CHECK(Constant::mn() == ParticleInfo(PID::neutron()).Mass());
    CHECK(Constant::mN() == (Constant::mp() + Constant::mn()) / 2.0);
    CHECK(Constant::mN2() == Constant::mN() * Constant::mN());
    CHECK(Constant::mlambda() == ParticleInfo(PID::lambda0()).Mass());
    CHECK(Constant::msigma0() == ParticleInfo(PID::sigma0()).Mass());
    // Sigma^- is 3112; PID::sigmam() is 3222, which the file lists as sigma+.
    CHECK(Constant::msigmam() == ParticleInfo(PID{3112}).Mass());

    // rho0 is in the database so nothing needs to hard-code it.
    CHECK(ParticleInfo(PID{113}).Name() == "rho0");
}
