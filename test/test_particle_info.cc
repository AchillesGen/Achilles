// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"

#include "Achilles/NuclearMass.hh"
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
        CHECK(info.Mass() == 10);
        CHECK(info.Width() == 1);
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

TEST_CASE("Nuclear PIDs are synthesized on demand", "[ParticleInfo]") {
    SECTION("Remnants absent from the database are built from the mass formula") {
        const auto pid = achilles::PID::MakeNucleus(6, 11);

        achilles::ParticleInfo info(pid);
        CHECK(info.ID() == pid);
        CHECK(info.IntCharge() == 18);
        CHECK(info.Name() == "C11");
        CHECK(info.Mass() == achilles::NuclearMass(6, 11));
        CHECK(info.IsNucleus());

        // The synthesized entry is cached, so a second lookup adds nothing.
        CHECK(achilles::ParticleInfo::Database().count(pid) == 1);
        const auto size = achilles::ParticleInfo::Database().size();
        achilles::ParticleInfo again(pid);
        CHECK(achilles::ParticleInfo::Database().size() == size);
        CHECK(again.Mass() == info.Mass());
    }

    SECTION("Database entries take precedence over the mass formula") {
        achilles::ParticleInfo carbon(achilles::PID::carbon());
        CHECK(carbon.Mass() == 11188);
        CHECK(carbon.Mass() != achilles::NuclearMass(6, 12));
        CHECK(carbon.IntCharge() == 18);
    }

    SECTION("Only well formed nuclear codes are accepted") {
        CHECK(achilles::PID::MakeNucleus(6, 12).valid_nucleus());
        CHECK(achilles::PID::MakeNucleus(1, 1).valid_nucleus());
        // More protons than nucleons
        CHECK_FALSE(achilles::PID::MakeNucleus(12, 6).valid_nucleus());
        // Not a nuclear code at all
        CHECK_FALSE(achilles::PID::proton().valid_nucleus());
        CHECK_THROWS_WITH(achilles::ParticleInfo(achilles::PID::MakeNucleus(12, 6)),
                          "Invalid PID: id=1000120060");
    }

    SECTION("Nuclear codes decode back to Z and A") {
        const auto pid = achilles::PID::MakeNucleus(18, 40);
        CHECK(pid.NuclearZ() == 18);
        CHECK(pid.NuclearA() == 40);
        CHECK(pid.NuclearL() == 0);
        CHECK(pid.NuclearI() == 0);
        CHECK(pid == achilles::PID::argon());
    }
}
