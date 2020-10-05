#include "nuchic/FInteractions.hh"
#include "nuchic/Particle.hh"

#include "fmt/format.h"
#include <iostream>

int main() {

    const std::string name = "Constant";
    nuchic::FortranInteraction interaction(name.c_str());
    nuchic::Particle part1{};
    nuchic::Particle part2{};

    double result = interaction.CrossSection(part1, part2);
    fmt::print("{}\n", result);

    std::array<double, 2> rans{0.5, 0.5};
    nuchic::ThreeVector mom = interaction.MakeMomentum(true, 100, rans); 
    fmt::print("{}\n", mom);
}
