#include "nuchic/Vegas.hh"
#include "nuchic/HardScattering.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Beams.hh"
#include "nuchic/Particle.hh"

#include "yaml-cpp/yaml.h"


class DummyDensity : public nuchic::Density {
    public:
        DummyDensity() = default;
        std::vector<nuchic::Particle> GetConfiguration() override {
            constexpr size_t Z = 6;
            nuchic::Particles particles;
            for(size_t i = 0; i < Z; ++i) {
                particles.emplace_back(nuchic::Particle(nuchic::PID::proton())); 
                particles.emplace_back(nuchic::Particle(nuchic::PID::neutron())); 
            } 
            return particles;
        }
};

nuchic::Histogram hist{50, 0, 400, "EnergyTransfer"};
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
   iterations: 30
   evaluations: 100000
    )node");

    YAML::Node config = YAML::Load(R"node(
QESettings:
    SpectralP: pke12_tot.data
    SpectralN: pke12_tot.data
    FermiGas: 1
    iform: 2
    )node");

    auto beam = std::make_shared<nuchic::Beam>(beams["Beams"].as<nuchic::Beam>());
    auto nucleus = std::make_shared<nuchic::Nucleus>(
            nuchic::Nucleus::MakeNucleus("12C", 0, 225, "c12.prova.txt",
                                         nuchic::Nucleus::FermiGasType::Global,
                                         std::make_unique<DummyDensity>()));

    auto mode = nuchic::RunMode::FixedAngle;
    //auto mode = nuchic::HardScatteringMode::FullPhaseSpace;

    //nuchic::FQESpectral hardScattering(config["QESettings"], beam, nucleus, mode);
    nuchic::FQEGlobalFermiGas hardScattering(config["QESettings"], beam, nucleus, mode);

    hardScattering.SetScatteringAngle(M_PI/180*15.0);

    nuchic::AdaptiveMap map(static_cast<size_t>(hardScattering.NVariables()));
    nuchic::Vegas vegas(map, node["Vegas"]);
    static constexpr double conv = 1e6; //conversion factor to obtain: nb/[MeV sr]
    auto xsec = [&](const std::vector<double> &x, const double &wgt) {
        // Generate the initial state nucleus
        nucleus -> GenerateConfig(); 
        auto particles = nucleus -> Nucleons();

        auto pswgt = hardScattering.GeneratePhaseSpace(particles, x);
        if(pswgt == 0) return pswgt;
        auto xsecwgt = hardScattering.CrossSection(particles);
        // fmt::print("theta = {:.5e}\tpswgt = {:.5e}\txsecwgt = {:.5e}\n",
        //            std::acos(cosTheta)*180/M_PI, pswgt, xsecwgt);
	if(fillHist) {
            auto omega = (particles[0].Momentum() - particles[1].Momentum()).E();
            auto niterations = node2["Vegas"]["iterations"].as<double>();
            hist.Fill(omega, pswgt*xsecwgt*wgt/niterations/(2*M_PI)*conv);
        }
        return pswgt*xsecwgt*conv; 
    };

    // vegas(xsec);
    // vegas.Clear();
    // vegas.Set(node2["Vegas"]);
    // hardScattering.SetHist(true);
    // fillHist = true;
    // vegas(xsec); 
    // hardScattering.GetHist().Save("test");
    // hist.Save("domega15_v3");
    // std::cout << hist.Integral() << std::endl;
}
