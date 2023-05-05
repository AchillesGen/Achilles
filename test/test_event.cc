#include "catch2/catch.hpp"
#include "catch_utils.hh"
#include "mock_classes.hh"

#include "Achilles/Beams.hh"
#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ProcessInfo.hh"

TEST_CASE("Initialize Event Parameters", "[Event]") {
    auto nuc = std::make_shared<MockNucleus>();
    auto beam = std::make_shared<MockBeam>();
    std::vector<double> rans{0};
    static constexpr achilles::FourVector lepton0{1000, 0, 0, 1000};
    static constexpr achilles::FourVector lepton1{313.073, 105.356, 174.207, -237.838};
    static constexpr achilles::FourVector hadron0{65.4247, 26.8702, -30.5306, -10.9449};
    static constexpr achilles::FourVector hadron1{1560.42, -78.4858, -204.738, 1226.89};

    std::vector<achilles::FourVector> moms = {lepton0, hadron0, lepton1, hadron1};
    achilles::Particles particles = {{achilles::PID::proton(), hadron0}};

    REQUIRE_CALL(*nuc, GenerateConfig()).TIMES(1);
    REQUIRE_CALL(*nuc, NNucleons()).LR_RETURN((12UL)).TIMES(1);
    static constexpr double vegas_wgt = 10;
    achilles::Event event(nuc, moms, vegas_wgt);

    SECTION("Nucleus and Beam set correctly") {
        CHECK(event.Momentum()[0] == lepton0);
        CHECK(event.Momentum()[1] == hadron0);
        CHECK(event.Momentum()[2] == lepton1);
        CHECK(event.Momentum()[3] == hadron1);
    }

    SECTION("Event can be finalized") {
        // Dummy carbon event
        achilles::Particles final = {
            {achilles::PID::proton(), hadron0, {}, achilles::ParticleStatus::initial_state},
            {achilles::PID::proton(), hadron1, {}, achilles::ParticleStatus::final_state},
            {achilles::PID::proton(), hadron0},
            {achilles::PID::proton(), hadron0},
            {achilles::PID::proton(), hadron0},
            {achilles::PID::proton(), hadron0},
            {achilles::PID::proton(), hadron0},
            {achilles::PID::neutron(), hadron0},
            {achilles::PID::neutron(), hadron0},
            {achilles::PID::neutron(), hadron0},
            {achilles::PID::neutron(), hadron0},
            {achilles::PID::neutron(), hadron0},
            {achilles::PID::neutron(), hadron0}};
        REQUIRE_CALL(*nuc, Nucleons()).LR_RETURN((final)).TIMES(AT_LEAST(13));

        event.Finalize();
        CHECK(event.Remnant().PID() == 1000050110);
        CHECK(event.Remnant().Mass() == 11 * achilles::Constant::mN);
    }
}
