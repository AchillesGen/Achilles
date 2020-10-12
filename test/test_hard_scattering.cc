#include "nuchic/Vegas.hh"
#include "nuchic/HardScattering.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Beams.hh"
#include "nuchic/Particle.hh"

#include "yaml-cpp/yaml.h"

constexpr size_t Z = 6;
nuchic::Particles Density() {
    nuchic::Particles particles;
    for(size_t i = 0; i < Z; ++i) {
        particles.emplace_back(nuchic::Particle(nuchic::PID::proton())); 
        particles.emplace_back(nuchic::Particle(nuchic::PID::neutron())); 
    } 
    return particles;
}

int main() {
    YAML::Node beams = YAML::Load(R"beam(
Beams:
    - Beam:
        - PID: 11
        - Beam Params:
            Type: Monochromatic
            Energy: 1108

    )beam"); 

    YAML::Node node = YAML::Load(R"node(
Vegas:
   seed: 123456789
   iterations: 10
   evaluations: 100000
    )node");

    YAML::Node node2 = YAML::Load(R"node(
Vegas:
   iterations: 20
   evaluations: 10000000
    )node");
    
    auto beam = beams["Beams"].as<nuchic::Beam>();
    auto nucleus = std::make_shared<nuchic::Nucleus>(
            nuchic::Nucleus::MakeNucleus("12C", 0, 225, "c12.prova.txt",
                                         nuchic::Nucleus::FermiGasType::Global, Density));

    nuchic::QESpectral hardScattering(beam, nucleus);

    nuchic::AdaptiveMap map(static_cast<size_t>(hardScattering.NVariables()));
    nuchic::Vegas vegas(map, node["Vegas"]);
    auto xsec = [&](const std::vector<double> &x, const double &wgt) {
        return hardScattering.Test(x, wgt); 
    };
    vegas(xsec);
    vegas.Clear();
    vegas.Set(node2["Vegas"]);
    hardScattering.SetHist(true);
    vegas(xsec);

    hardScattering.GetHist().Save("test");
}
