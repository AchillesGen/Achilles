#include "catch2/catch.hpp"

#include "nuchic/ParticleInfo.hh"

TEST_CASE("ParticleInfo", "[ParticleInfo]") {
    SECTION("Must be a valid particle") {
        CHECK_THROWS_WITH(nuchic::ParticleInfo(23413), 
                          "Invalid PID: id=23413");

        CHECK_THROWS_WITH(nuchic::ParticleInfo(nuchic::PID(23413)),
                          "Invalid PID: id=23413");
       
        for(const auto &entry : nuchic::ParticleInfo::Database()) {
            CHECK_NOTHROW(nuchic::ParticleInfo(entry.first));
        }
    }

    SECTION("Add Entry to database") {
        const auto size = nuchic::ParticleInfo::Database().size();
        auto entry = std::make_shared<nuchic::ParticleInfoEntry>(nuchic::PID(123456789), 0, 0, 0, 0,
                                                                 0, 0, 0, true, false, "test", "anti-test");
        nuchic::ParticleInfo info(entry);

        CHECK(nuchic::ParticleInfo::Database().size() == size + 1);
    }

    SECTION("Negative PIDs are interpreted as anti-particles") {
        nuchic::ParticleInfo info1(nuchic::PID::proton()), info2(-2212); 

        CHECK(info1.Anti() == info2);
        CHECK(info1.Anti().Anti() == info1);
    }

    SECTION("Valid properties are returned") {
        // dummy particle
        auto entry = std::make_shared<nuchic::ParticleInfoEntry>(nuchic::PID(123456789), 10, 1, 1, 1,
                                                                 1, 1, 0, true, false, "test", "anti-test");
        nuchic::ParticleInfo info(entry);
        nuchic::ParticleInfo ainfo(entry, true);

        CHECK(info.Name() == "test");
        CHECK(ainfo.Name() == "anti-test");

        CHECK(info.ID() == nuchic::PID(123456789));
        CHECK(ainfo.ID() == nuchic::PID(123456789));

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
        CHECK(info.Charge() == 1./3);
        CHECK(ainfo.Charge() == -1./3);
        CHECK(info.IntSpin() == 1);
        CHECK(info.Spin() == 1./2);

        CHECK(info.SelfAnti() == false);
        CHECK(info.Majorana() == false);
        CHECK(info.Stable() == 1);
        CHECK(info.IsStable() == true);
        CHECK(info.IsMassive() == true);
        CHECK(info.Mass() == 10);
        CHECK(info.Width() == 1);
    }

    SECTION("Equality of two particles") {
        nuchic::ParticleInfo info1(nuchic::PID::proton()), info2(nuchic::PID::proton()),
                             info3(nuchic::PID::neutron()), info4(-2212);

        // Only same PIDs are equal
        CHECK(info1 == info2);
        CHECK(info1 != info3);
        // Anti-particles are not equal to particles
        CHECK(info1 != info4);
    }
}
