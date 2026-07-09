#include "catch2/catch.hpp"

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

TEST_CASE("Particle Containers", "[ParticleInfo]") {
    achilles::ParticleInfo::RegisterContainer(achilles::PID::nucleon(), "nucleon", {2212, 2112});
    achilles::ParticleInfo info(achilles::PID::nucleon());
    achilles::ParticleInfo proton(achilles::PID::proton());
    achilles::ParticleInfo neutron(achilles::PID::neutron());

    SECTION("Creates groups") {
        CHECK(info.IsGroup());
    }

    SECTION("Group Membership") {
        CHECK(info.Size() == 2);
        CHECK(info[0] == proton);
        CHECK(info[1] == neutron);
        CHECK(info.Includes(achilles::PID::proton()));
        CHECK(info.Includes(achilles::PID::neutron()));
    }

    SECTION("Decompose") {
        auto ids = info.Decompose();
        CHECK(ids.size() == 2);
        CHECK(ids[0] == proton);
        CHECK(ids[1] == neutron);
    }

    SECTION("Properties") {
        CHECK(info.Mass() == proton.Mass());
    }

    SECTION("Masses must be the same") {
        CHECK_THROWS_AS(achilles::ParticleInfo::RegisterContainer(achilles::PID::lepton(), "lepton",
                                                                  {11, 13, 15}),
                        std::runtime_error);
    }
}

TEST_CASE("FinalizeContainers", "[ParticleInfo]") {
    using achilles::ParticleInfo;
    using achilles::PID;

    // Massless (light) quarks only, excluding the undefined placeholder (PID 0).
    auto light_quarks = [](const ParticleInfo &p) {
        return p.ID() != PID::undefined() && p.IsQuark() && !p.IsMassive();
    };

    SECTION("Rule-based population from the table") {
        ParticleInfo::RegisterContainerRule(PID::quark(), "quark", light_quarks);

        // Empty until finalized.
        CHECK(ParticleInfo(PID::quark()).Size() == 0);

        ParticleInfo::FinalizeContainers();

        ParticleInfo quark(PID::quark());
        CHECK(quark.IsGroup());
        CHECK(quark.Size() == 3);
        CHECK(quark.Includes(PID::down()));
        CHECK(quark.Includes(PID::up()));
        CHECK(quark.Includes(PID::strange()));
        CHECK_FALSE(quark.Includes(PID::charm())); // massive -> excluded
        CHECK_FALSE(quark.Includes(PID::gluon())); // not a quark
        CHECK(quark.Mass() == 0);                  // adopted massless
    }

    SECTION("Idempotent") {
        ParticleInfo::RegisterContainerRule(PID::quark(), "quark", light_quarks);
        ParticleInfo::FinalizeContainers();
        ParticleInfo::FinalizeContainers();
        CHECK(ParticleInfo(PID::quark()).Size() == 3); // not doubled
    }

    SECTION("BSM/open-world: re-finalizing tracks the table") {
        ParticleInfo::RegisterContainerRule(PID::quark(), "quark", light_quarks);
        ParticleInfo::FinalizeContainers();
        REQUIRE(ParticleInfo(PID::quark()).Size() == 3);

        // A newly registered massless quark-like state is picked up on re-finalize.
        const PID bsm_quark{7};
        ParticleInfo::Database()[bsm_quark] = std::make_shared<achilles::ParticleInfoEntry>(
            bsm_quark, 0, 0, -1, 3, 1, 1, 2, 0, false, false, "d'", "d'bar");
        ParticleInfo::FinalizeContainers();
        CHECK(ParticleInfo(PID::quark()).Size() == 4);
        CHECK(ParticleInfo(PID::quark()).Includes(bsm_quark));

        // Clean up so the added state does not leak into other test cases.
        ParticleInfo::Database().erase(bsm_quark);
        ParticleInfo::FinalizeContainers();
        CHECK(ParticleInfo(PID::quark()).Size() == 3);
    }

    SECTION("Containers are never nested") {
        // A rule that would match container PIDs (90-95) must still yield an empty
        // container: FinalizeContainers skips groups and rule containers (incl. self).
        ParticleInfo::RegisterContainerRule(PID::quark(), "quark", [](const ParticleInfo &p) {
            return p.ID().AsInt() >= 90 && p.ID().AsInt() <= 95;
        });
        ParticleInfo::FinalizeContainers();
        CHECK(ParticleInfo(PID::quark()).Size() == 0);
    }
}
