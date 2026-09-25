// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_test_macros.hpp"

#include "Achilles/Configuration.hh"
#include "Achilles/Particle.hh"
#include "Achilles/PhysicalUnits.hh"

TEST_CASE("DensityConfiguration", "[Configuration]") {
    achilles::DensityConfiguration config("data/configurations/QMC_configs.out.gz");
    auto particles = config.GetConfiguration();
    CHECK(particles.size() == 12);
    size_t nproton = 0, nneutron = 0;
    for(const auto &particle : particles) {
        if(particle.ID() == achilles::PID::proton()) nproton++;
        if(particle.ID() == achilles::PID::neutron()) nneutron++;
    }
    CHECK(nproton == 6);
    CHECK(nneutron == 6);

    // The file gives the coordinates in fm: a 12C nucleon sits a few fm out,
    // so reading them in any other unit is off by hbar*c and shows up here.
    bool any_off_origin = false;
    for(const auto &particle : particles) {
        const double r = particle.Position().Magnitude().in(achilles::units::fm);
        CHECK(r < 15);
        if(r > 0.1) any_off_origin = true;
    }
    CHECK(any_off_origin);
}
