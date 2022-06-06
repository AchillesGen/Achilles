#include "Achilles/FInteractions.hh"
#include "Achilles/Particle.hh"

#include "fmt/format.h"
#include <iostream>

int main() {

    const YAML::Node name = YAML::Load("Name: Constant");
    achilles::FortranInteraction interaction(name);
    achilles::Particle part1{};
    achilles::Particle part2{};

    double result = interaction.CrossSection(part1, part2);
    fmt::print("{}\n", result);

    std::array<double, 2> rans{0.5, 0.5};
    achilles::ThreeVector mom = interaction.MakeMomentum(true, 100, rans); 
    fmt::print("{}\n", mom);
}
