#include "catch2/catch.hpp"
#include "mock_classes.hh"

#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ProcessInfo.hh"

using achilles::Particle;

TEST_CASE("Initialize Event Parameters", "[Event]") {
    auto nuc = std::make_shared<MockNucleus>();
    auto beam = std::make_shared<MockBeam>();
    std::vector<double> rans{0};
    static constexpr achilles::FourVector lepton0{1000, 0, 0, 1000};
    static constexpr achilles::FourVector lepton1{313.073, 105.356, 174.207, -237.838};
    static constexpr achilles::FourVector hadron0{65.4247, 26.8702, -30.5306, -10.9449};
    static constexpr achilles::FourVector hadron1{1560.42, -78.4858, -204.738, 1226.89};

    std::vector<achilles::FourVector> moms = {lepton0, hadron0, lepton1, hadron1};
    achilles::Particles particles = {Particle{achilles::PID::proton(), hadron0}};

    REQUIRE_CALL(*nuc, GenerateConfig()).TIMES(1).RETURN((particles));
    static constexpr double vegas_wgt = 10;
    achilles::Event event(nuc, moms, vegas_wgt);

    SECTION("Nucleus and Beam set correctly") {
        CHECK(event.Momentum()[0] == lepton0);
        CHECK(event.Momentum()[1] == hadron0);
        CHECK(event.Momentum()[2] == lepton1);
        CHECK(event.Momentum()[3] == hadron1);
    }
}

TEST_CASE("Finalize Event", "[Event]") {
    static constexpr achilles::FourVector hadron0{65.4247, 26.8702, -30.5306, -10.9449};
    static constexpr achilles::FourVector hadron1{1560.42, -78.4858, -204.738, 1226.89};

    // Dummy carbon event
    achilles::Particles final = {
        Particle{achilles::PID::proton(), hadron0, {}, achilles::ParticleStatus::initial_state},
        Particle{achilles::PID::proton(), hadron1, {}, achilles::ParticleStatus::final_state},
        Particle{achilles::PID::proton(), hadron0},
        Particle{achilles::PID::proton(), hadron0},
        Particle{achilles::PID::proton(), hadron0},
        Particle{achilles::PID::proton(), hadron0},
        Particle{achilles::PID::proton(), hadron0},
        Particle{achilles::PID::neutron(), hadron0},
        Particle{achilles::PID::neutron(), hadron0},
        Particle{achilles::PID::neutron(), hadron0},
        Particle{achilles::PID::neutron(), hadron0},
        Particle{achilles::PID::neutron(), hadron0},
        Particle{achilles::PID::neutron(), hadron0}};

    achilles::Event event{};
    event.Hadrons() = final;
    event.Finalize();
    CHECK(event.Remnant().PID() == 1000050110);
    CHECK(event.Remnant().Mass() == 11 * achilles::Constant::mN);
}
