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

nuchic::Histogram hist{200, 0, 1000, "EnergyTransfer"};
bool fillHist{false};

int main() {
    YAML::Node beams = YAML::Load(R"beam(
Beams:
    - Beam:
        - PID: 11
        - Beam Params:
            Type: Monochromatic
            Energy: 1000

    )beam"); 

    YAML::Node node = YAML::Load(R"node(
Vegas:
   seed: 123456789
   iterations: 20
   evaluations: 100000
    )node");

    YAML::Node node2 = YAML::Load(R"node(
Vegas:
   iterations: 20
   evaluations: 10000000
    )node");

    YAML::Node config = YAML::Load(R"node(
QESettings:
    SpectralP: pke12_tot.data
    SpectralN: pke12_tot.data
    FermiGas: 0
    iform: 2
    )node");

    auto beam = beams["Beams"].as<nuchic::Beam>();
    auto nucleus = std::make_shared<nuchic::Nucleus>(
            nuchic::Nucleus::MakeNucleus("12C", 0, 225, "c12.prova.txt",
                                         nuchic::Nucleus::FermiGasType::Global, Density));

    auto mode = nuchic::HardScatteringMode::FixedAngle;
    nuchic::FQESpectral hardScattering(config["QESettings"], beam, nucleus, mode);
    hardScattering.SetScatteringAngle(M_PI/180*15.0);

    nuchic::AdaptiveMap map(static_cast<size_t>(hardScattering.NVariables()));
    nuchic::Vegas vegas(map, node["Vegas"]);
    auto xsec = [&](const std::vector<double> &x, const double &wgt) {
        auto particles = hardScattering.GeneratePhaseSpace(x);
        if(particles[2].E() == 0) return 0.0;
        double cosTheta = particles[1].Momentum().CosAngle(particles[0].Momentum());
        if(std::abs(cosTheta) > std::cos(M_PI/180)) return 0.0;
        double pswgt = hardScattering.PhaseSpaceWeight(particles);
        if(pswgt == 0) return pswgt;
        double xsecwgt = hardScattering.CrossSection(particles);
        // fmt::print("theta = {:.5e}\tpswgt = {:.5e}\txsecwgt = {:.5e}\n",
        //            std::acos(cosTheta)*180/M_PI, pswgt, xsecwgt);
        if(fillHist) {
            auto omega = (particles[0].Momentum() - particles[1].Momentum()).E();
            hist.Fill(omega, pswgt*xsecwgt*wgt);
        }
        return pswgt*xsecwgt; 
    };
    vegas(xsec);
    vegas.Clear();
    vegas.Set(node2["Vegas"]);
    hardScattering.SetHist(true);
    fillHist = true;
    vegas(xsec);

    hardScattering.GetHist().Save("test");
    hist.Save("cosTheta");
}
