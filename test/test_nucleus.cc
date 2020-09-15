#include <iostream>

#include "catch2/catch.hpp"
#include "catch2/trompeloeil.hpp"

#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"

using Particles = std::vector<nuchic::Particle>;

Particles density(size_t Z, size_t A) {
    Particles parts;
    for(size_t i = 0; i < A; ++i) {
        nuchic::PID pid = i < Z ? nuchic::PID::proton() : nuchic::PID::neutron();
        parts.push_back(nuchic::Particle(pid));
    }
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
    const std::string radius = "c12.prova.txt";

    SECTION("Nucleus must have more nucleons than protons") {
        auto gasType = nuchic::Nucleus::FermiGasType::Global;

        static constexpr std::size_t Z = 10, A = 15;
        auto fdensity = [&](){ return density(Z, A); };
        nuchic::Nucleus nuc(Z, A, 0, 0, radius, gasType, fdensity);

        CHECK(nuc.NNucleons() == A);
        CHECK(nuc.NProtons() == Z);
        CHECK(nuc.NNeutrons() == A-Z);

        CHECK_THROWS(nuchic::Nucleus(A, Z, 0, 0, radius, gasType, fdensity));

        CHECK(nuc.Radius() > 0);
        CHECK(nuc.PotentialEnergy() >= 0);
    }

    SECTION("MakeNucleus creates a proper Nucleus", "[nucleus]") {
        auto gasType = nuchic::Nucleus::FermiGasType::Global;

        for(auto name : names) {
            auto fdensity = [&](){ return density(name.second, 2*name.second); };
            auto nuc = nuchic::Nucleus::MakeNucleus(name.first, 0, 0, radius, gasType, fdensity);
            CHECK(nuc.NNucleons() == name.second*2);
            CHECK(nuc.NProtons() == name.second);
            CHECK(nuc.NNeutrons() == name.second);
        }
    }
}

static constexpr double maxRadius = 5.0;
static constexpr size_t nTrials = 50;

TEST_CASE("Nucleus has proper Fermi momentum", "[nucleus]") {
    const std::string radiusFile = "c12.prova.txt";
    constexpr double kf = 225;

    SECTION("Fermi momentum is correct for global Fermi gas") {
        auto gasType = nuchic::Nucleus::FermiGasType::Global;

        static constexpr std::size_t Z = 6, A = 12;
        auto fdensity = [&]{ return density(Z, A); };
        nuchic::Nucleus nuc(Z, A, 0, kf, radiusFile, gasType, fdensity);

        auto radius = GENERATE(take(nTrials, random(0.0, maxRadius)));
        CHECK(nuc.FermiMomentum(radius) == kf);
    }

    SECTION("Fermi momentum is correct for local Fermi gas") {
        auto gasType = nuchic::Nucleus::FermiGasType::Local;
        static constexpr double maxRho = 0.09;
        const double maxFermiMomentum = std::cbrt(maxRho*3*M_PI*M_PI)*nuchic::Constant::HBARC;

        static constexpr std::size_t Z = 6, A = 12;
        auto fdensity = [&]{ return density(Z, A); };
        nuchic::Nucleus nuc(Z, A, 0, kf, radiusFile, gasType, fdensity);

        auto radius = GENERATE(take(nTrials, random(0.0, maxRadius)));
        CHECK(nuc.FermiMomentum(radius) <= maxFermiMomentum);
    }
}

TEST_CASE("Nucleus generates proper momentum", "[nucleus]") {
    const std::string radiusFile = "c12.prova.txt";
    constexpr double kf = 225;

    SECTION("Global Fermi gas") {
        auto gasType = nuchic::Nucleus::FermiGasType::Global;

        static constexpr std::size_t Z = 6, A = 12;
        auto fdensity = [&]{ return density(Z, A); };
        nuchic::Nucleus nuc(Z, A, 0, kf, radiusFile, gasType, fdensity);

        auto radius = GENERATE(take(nTrials, random(0.0, maxRadius)));
        auto mom = nuc.GenerateMomentum(radius);
        CHECK(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2] <= kf*kf);
    }

    SECTION("Local Fermi gas") {
        auto gasType = nuchic::Nucleus::FermiGasType::Local;
        static constexpr double maxRho = 0.09;
        const double maxFermiMomentum = std::cbrt(maxRho*3*M_PI*M_PI)*nuchic::Constant::HBARC;

        static constexpr std::size_t Z = 6, A = 12;
        auto fdensity = [&]{ return density(Z, A); };
        nuchic::Nucleus nuc(Z, A, 0, kf, radiusFile, gasType, fdensity);

        auto radius = GENERATE(take(nTrials, random(0.0, maxRadius)));
        auto mom = nuc.GenerateMomentum(radius);
        CHECK(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2] <= pow(maxFermiMomentum, 2));
    }
}
