#include "catch2/catch.hpp"

#include "Achilles/Configuration.hh"
#include "Achilles/Particle.hh"

TEST_CASE("DensityConfiguration", "[Configuration]") {
    achilles::DensityConfiguration config("data/configurations/QMC_configs.out.gz"); 
    auto particles = config.GetConfiguration();
    CHECK(particles.size() == 12);
    size_t nproton=0, nneutron=0;
    for(const auto &particle : particles) {
        if(particle.ID() == achilles::PID::proton()) nproton++;
        if(particle.ID() == achilles::PID::neutron()) nneutron++;
    }
    CHECK(nproton == 6);
    CHECK(nneutron == 6);
}
