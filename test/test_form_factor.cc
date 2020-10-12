#include "catch2/catch.hpp"

#include <iostream>

#include "nuchic/FormFactor.hh"
#include "nuchic/Units.hh"

#include "yaml-cpp/yaml.h"

using nuchic::operator""_GeV;

TEST_CASE("Dipole", "[FormFactor]") {
    YAML::Node formFactorYAML = YAML::LoadFile("FormFactors.yml");
    nuchic::Dipole dipole(formFactorYAML["Dipole"]);
    nuchic::Kelly kelly(formFactorYAML["Kelly"]);
    nuchic::BBBA bbba(formFactorYAML["BBBA"]);
    nuchic::ArringtonHill hill(formFactorYAML["ArringtonHill"]);

    double Q2{};
    for(size_t i = 1; i < 100; ++i) {
        Q2 = static_cast<double>(i*i*100);
        std::cout << sqrt(Q2) << " " << dipole(Q2/1_GeV/1_GeV).Gep;
        std::cout << " " << kelly(Q2/1_GeV/1_GeV).Gep;
        std::cout << " " << bbba(Q2/1_GeV/1_GeV).Gep;
        std::cout << " " << hill(Q2/1_GeV/1_GeV).Gep << std::endl;
    }
}
