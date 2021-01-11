#include <iostream>

#include "catch2/catch.hpp"
#include "catch2/trompeloeil.hpp"

#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"

#include "catch_utils.hh"

std::string const& RandomNucleusGenerator::get() const {
    return current_string;
}

// This helper function provides a nicer UX when instantiating the generator
// Notice that it returns an instance of GeneratorWrapper<std::string>, which
// is a value-wrapper around std::unique_ptr<IGenerator<std::string>>.
Catch::Generators::GeneratorWrapper<std::string> randomNucleus(size_t length) {
    return Catch::Generators::GeneratorWrapper<std::string>(std::unique_ptr<Catch::Generators::IGenerator<std::string>>(new RandomNucleusGenerator(length)));
}

class MockDensity : public trompeloeil::mock_interface<nuchic::Density> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(GetConfiguration);
};

const std::string dFile = "data/c12.prova.txt";

TEST_CASE("Nucleus construction", "[nucleus]") {
    auto fermiGas = GENERATE(nuchic::Nucleus::FermiGasType::Local,
                             nuchic::Nucleus::FermiGasType::Global);

    SECTION("Nucleus must have more nucleons than protons") {
        static constexpr std::size_t Z = 6, A = 12;
        nuchic::Particles particles;
        for(size_t i = 0; i < Z; ++i) {
            particles.emplace_back(nuchic::PID::proton());
            particles.emplace_back(nuchic::PID::neutron());
        }

        auto density1 = std::make_unique<MockDensity>();
        REQUIRE_CALL(*density1, GetConfiguration())
            .TIMES(1)
            .RETURN(particles);
        CHECK_NOTHROW(nuchic::Nucleus(Z, A, 0, 0, dFile, fermiGas, std::move(density1)));

        auto density2 = std::make_unique<MockDensity>();
        REQUIRE_CALL(*density2, GetConfiguration())
            .TIMES(1)
            .RETURN(particles);
        nuchic::Nucleus nuc(Z, A, 0, 0, dFile, fermiGas, std::move(density2));

        CHECK(nuc.NNucleons() == A);
        CHECK(nuc.NProtons() == Z);
        CHECK(nuc.NNeutrons() == A-Z);
        CHECK(nuc.Radius() > 0);
        // CHECK(nuc.PotentialEnergy() > 0);

        auto density3 = std::make_unique<MockDensity>();
        REQUIRE_CALL(*density3, GetConfiguration())
            .TIMES(0);
        std::string errorMsg = "Requires the number of protons to be less than the total";
        errorMsg += " number of nucleons. Got " + std::to_string(A);
        errorMsg += " protons and " + std::to_string(Z) + " nucleons";
        CHECK_THROWS_WITH(nuchic::Nucleus(A, Z, 0, 0, dFile, fermiGas, std::move(density3)),
                          errorMsg);
    }

    SECTION("Nucleus needs valid density file") {
        auto density = std::make_unique<MockDensity>();
        REQUIRE_CALL(*density, GetConfiguration())
            .TIMES(0);
        static constexpr std::size_t Z = 6, A = 12;

        CHECK_THROWS_WITH(nuchic::Nucleus(Z, A, 0, 0, "dummy.txt", fermiGas, std::move(density)),
                          "Nucleus: Density file dummy.txt does not exist.");
    }

    SECTION("Density must produce correct number of protons and neutrons") {
        static constexpr std::size_t Z = 6, A = 12;
        nuchic::Particles particles;
        for(size_t i = 0; i < Z; ++i) {
            particles.emplace_back(nuchic::PID::proton());
            particles.emplace_back(nuchic::PID::neutron());
        }

        auto density1 = std::make_unique<MockDensity>();
        REQUIRE_CALL(*density1, GetConfiguration())
            .TIMES(1)
            .RETURN(particles);
        CHECK_NOTHROW(nuchic::Nucleus(Z, A, 0, 0, dFile, fermiGas, std::move(density1)));

        auto density2 = std::make_unique<MockDensity>();
        REQUIRE_CALL(*density2, GetConfiguration())
            .TIMES(1)
            .RETURN(particles);
        CHECK_THROWS_WITH(nuchic::Nucleus(Z, A+1, 0, 0, dFile, fermiGas, std::move(density2)),
                          "Invalid density function! Incorrect number of nucleons.");

        auto density3 = std::make_unique<MockDensity>();
        REQUIRE_CALL(*density3, GetConfiguration())
            .TIMES(1)
            .RETURN(particles);
        CHECK_THROWS_WITH(nuchic::Nucleus(Z+1, A, 0, 0, dFile, fermiGas, std::move(density3)),
                          "Invalid density function! Incorrect number of protons or neutrons.");
    }
}

TEST_CASE("Make Nucleus", "[Nucleus]") {
    const auto fermiGas = nuchic::Nucleus::FermiGasType::Local;

    SECTION("Creates a proper Nucleus") {
        auto name = GENERATE(table<std::string, size_t>({
                    {"2H",    1},
                    {"4He",   2},
                    {"6Li",   3},
                    {"12C",   6},
                    {"16O",   8}, 
                    {"26Al", 13}, 
                    {"36Ar", 18}, 
                    {"40Ca", 20},
                    {"52Fe", 26}
                }));

        nuchic::Particles particles;
        for(size_t i = 0; i < std::get<1>(name); ++i) {
            particles.emplace_back(nuchic::PID::proton());
            particles.emplace_back(nuchic::PID::neutron());
        }

        auto density = std::make_unique<MockDensity>();
        REQUIRE_CALL(*density, GetConfiguration())
            .TIMES(1)
            .RETURN(particles);

        auto nuc = nuchic::Nucleus::MakeNucleus(std::get<0>(name), 0, 0, dFile, fermiGas, std::move(density));
        CHECK(nuc.NNucleons() == 2*std::get<1>(name));
        CHECK(nuc.NProtons() == std::get<1>(name));
        CHECK(nuc.NNeutrons() == std::get<1>(name));
    }

    SECTION("Throws on invalid element") {
        auto name = GENERATE(take(30, randomNucleus(200)));
        auto density = std::make_unique<MockDensity>();
        FORBID_CALL(*density, GetConfiguration());

        spdlog::info("Nucleus name: {}", name);
        const std::regex regex("([0-9]+)([a-zA-Z]+)");
        std::smatch match;

        std::regex_match(name, match, regex);
        CHECK_THROWS_WITH(nuchic::Nucleus::MakeNucleus(name, 0, 0, dFile, fermiGas, std::move(density)),
                          fmt::format("Invalid nucleus: {} does not exist.", match[2]));
    }
}
