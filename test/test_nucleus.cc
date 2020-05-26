#include <iostream>

#include "catch2/catch.hpp"
#include "catch2/trompeloeil.hpp"

#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"

nuchic::Particles density() {
    nuchic::Particles parts;
    parts.push_back(nuchic::Particle());
    return parts;
}

const std::map<std::string, std::size_t> names = {
    {"2H",    1},
    {"4He",   2},
    {"6Li",   3},
    {"12C",   6},
    {"16O",   8}, 
    {"26Al", 13}, 
    {"36Ar", 18}, 
    {"40Ca", 20},
    {"52Fe", 26}
};

TEST_CASE("Nucleus is constructed properly", "[nucleus]") {
    SECTION("Nucleus must have more nucleons than protons") {
        std::size_t Z = 10, A = 15;
        nuchic::Nucleus nuc(Z, A, 0, 0, density);

        CHECK(nuc.NNucleons() == A);
        CHECK(nuc.NProtons() == Z);
        CHECK(nuc.NNeutrons() == A-Z);

        CHECK_THROWS(nuchic::Nucleus(A, Z, 0, 0, density));

        CHECK(nuc.Radius() > 0);
        CHECK(nuc.PotentialEnergy() > 0);
    }

    SECTION("MakeNucleus creates a proper Nucleus", "[nucleus]") {
        for(auto name : names) {
            auto nuc = nuchic::Nucleus::MakeNucleus(name.first, 0, 0, density);
            CHECK(nuc.NNucleons() == name.second*2);
            CHECK(nuc.NProtons() == name.second);
            CHECK(nuc.NNeutrons() == name.second);
        }
    }
}
