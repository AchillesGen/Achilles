#include "catch2/catch.hpp"
#include "mock_classes.hh"

#include <sstream>

#include "plugins/NuHepMC/NuHepMCWriter.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Version.hh"

TEST_CASE("Passes Validator", "[NuHepMC]") {
    // Setup dummy event
    static constexpr achilles::FourVector hadron0{65.4247, 26.8702, -30.5306, -10.9449};
    static constexpr achilles::FourVector hadron1{1560.42, -78.4858, -204.738, 1226.89};
    achilles::Particles hadrons = {
        {achilles::PID::proton(), hadron0, {}, achilles::ParticleStatus::initial_state},
        {achilles::PID::proton(), hadron1, {}, achilles::ParticleStatus::final_state}};
    achilles::Particles leptons = {
        {achilles::PID::electron(), hadron0, {}, achilles::ParticleStatus::initial_state},
        {achilles::PID::electron(), hadron1, {}, achilles::ParticleStatus::final_state}};
    achilles::NuclearRemnant remnant(11, 5);

    auto nuc = std::make_shared<MockNucleus>();
    auto id = achilles::PID::carbon();
    REQUIRE_CALL(*nuc, ID())
        .TIMES(AT_LEAST(1))
        .LR_RETURN(id);
    
    const MockEvent event;
    double wgt = 1.0;
    REQUIRE_CALL(event, Hadrons())
        .TIMES(1)
        .LR_RETURN((hadrons));
    REQUIRE_CALL(event, Leptons())
        .TIMES(1)
        .LR_RETURN((leptons));
    REQUIRE_CALL(event, CurrentNucleus())
        .TIMES(AT_LEAST(1))
        .LR_RETURN((nuc));
    REQUIRE_CALL(event, Remnant())
        .TIMES(AT_LEAST(1))
        .LR_RETURN((remnant));
    REQUIRE_CALL(event, Weight())
        .TIMES(AT_LEAST(1))
        .LR_RETURN((wgt));

    // Force saving of HepMC file
    {
        achilles::NuHepMCWriter writer("data.hepmc", false);
        writer.WriteHeader("dummy");

        writer.Write(event);
    }
    auto result = std::system("./bin/NuHepMCReferenceValidator data.hepmc");
    CHECK(result == 0);
    // result = std::system("rm data.hepmc");
}
