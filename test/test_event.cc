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

TEST_CASE("Test Event Copy Constructor", "[Event]") {
    auto nuc = std::make_shared<MockNucleus>();
    static constexpr achilles::FourVector lepton0{1000, 0, 0, 1000};
    static constexpr achilles::FourVector lepton1{313.073, 105.356, 174.207, -237.838};
    static constexpr achilles::FourVector hadron0{65.4247, 26.8702, -30.5306, -10.9449};
    static constexpr achilles::FourVector hadron1{1560.42, -78.4858, -204.738, 1226.89};

    std::vector<achilles::FourVector> moms = {lepton0, hadron0, lepton1, hadron1};
    achilles::Particles particles = {Particle{achilles::PID::proton(), hadron0}};

    REQUIRE_CALL(*nuc, GenerateConfig()).TIMES(1).RETURN(particles);
    static constexpr double vegas_wgt = 10;
    achilles::Event event(nuc, moms, vegas_wgt);

    // Copy constructor test
    achilles::Event copied_event = event;

    SECTION("Check if copy constructor works correctly") {
        CHECK(copied_event.Momentum()[0] == lepton0);
        CHECK(copied_event.Momentum()[1] == hadron0);
        CHECK(copied_event.Momentum()[2] == lepton1);
        CHECK(copied_event.Momentum()[3] == hadron1);
    }

    SECTION("Ensure original and copied events are independent") {
        // Modify the original event and ensure it does not affect the copied event
        event.Momentum()[0] = {0, 0, 0, 0}; // Modify first momentum of the original event

        CHECK(event.Momentum()[0] !=
              copied_event.Momentum()[0]); // The copied event should remain unchanged
        CHECK(event.Momentum()[1] ==
              copied_event.Momentum()[1]); // Other parts should remain the same
    }
}

TEST_CASE("Test Event Assignment Operator", "[Event]") {
    auto nuc = std::make_shared<MockNucleus>();
    static constexpr achilles::FourVector lepton0{1000, 0, 0, 1000};
    static constexpr achilles::FourVector lepton1{313.073, 105.356, 174.207, -237.838};
    static constexpr achilles::FourVector hadron0{65.4247, 26.8702, -30.5306, -10.9449};
    static constexpr achilles::FourVector hadron1{1560.42, -78.4858, -204.738, 1226.89};

    std::vector<achilles::FourVector> moms = {lepton0, hadron0, lepton1, hadron1};
    achilles::Particles particles = {Particle{achilles::PID::proton(), hadron0}};

    REQUIRE_CALL(*nuc, GenerateConfig()).TIMES(1).RETURN(particles);
    static constexpr double vegas_wgt = 10;
    achilles::Event event(nuc, moms, vegas_wgt);

    // Assignment operator test
    achilles::Event assigned_event;
    assigned_event = event;

    SECTION("Check if assignment operator works correctly") {
        CHECK(assigned_event.Momentum()[0] == lepton0);
        CHECK(assigned_event.Momentum()[1] == hadron0);
        CHECK(assigned_event.Momentum()[2] == lepton1);
        CHECK(assigned_event.Momentum()[3] == hadron1);
    }

    SECTION("Ensure original and assigned events are independent") {
        // Modify the original event and ensure it does not affect the assigned event
        event.Momentum()[0] = {0, 0, 0, 0}; // Modify first momentum of the original event

        CHECK(event.Momentum()[0] !=
              assigned_event.Momentum()[0]); // The assigned event should remain unchanged
        CHECK(event.Momentum()[1] ==
              assigned_event.Momentum()[1]); // Other parts should remain the same
    }
}
