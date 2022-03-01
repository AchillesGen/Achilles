#include "catch2/catch.hpp"

#include "nuchic/Configuration.hh"
#include "nuchic/Particle.hh"

TEST_CASE("DensityConfiguration", "[Configuration]") {
    nuchic::DensityConfiguration config("data/configurations/QMC_configs.out.gz"); 
    auto particles = config.GetConfiguration();
    CHECK(particles.size() == 12);
    size_t nproton=0, nneutron=0;
    for(const auto &particle : particles) {
        if(particle.ID() == nuchic::PID::proton()) nproton++;
        if(particle.ID() == nuchic::PID::neutron()) nneutron++;
    }
    CHECK(nproton == 6);
    CHECK(nneutron == 6);
}
