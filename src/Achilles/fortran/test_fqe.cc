#include "Achilles/Vegas.hh"
#include "Achilles/HardScattering.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Particle.hh"

#include "yaml-cpp/yaml.h"


class DummyDensity : public achilles::Density {
    public:
        DummyDensity() = default;
        std::vector<achilles::Particle> GetConfiguration() override {
            constexpr size_t Z = 6;
            achilles::Particles particles;
            for(size_t i = 0; i < Z; ++i) {
                particles.emplace_back(achilles::Particle(achilles::PID::proton())); 
                particles.emplace_back(achilles::Particle(achilles::PID::neutron())); 
            } 
            return particles;
        }
};

achilles::Histogram hist{50, 0, 400, "EnergyTransfer"};
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

    auto beam = std::make_shared<achilles::Beam>(beams["Beams"].as<achilles::Beam>());
    auto nucleus = std::make_shared<achilles::Nucleus>(
            achilles::Nucleus::MakeNucleus("12C", 0, 225, "c12.prova.txt",
                                         achilles::Nucleus::FermiGasType::Global,
                                         std::make_unique<DummyDensity>()));

    auto mode = achilles::RunMode::FixedAngle;
    //auto mode = achilles::HardScatteringMode::FullPhaseSpace;

    //achilles::FQESpectral hardScattering(config["QESettings"], beam, nucleus, mode);
    achilles::FQEGlobalFermiGas hardScattering(config["QESettings"], beam, nucleus, mode);

    hardScattering.SetScatteringAngle(M_PI/180*15.0);

    achilles::AdaptiveMap map(static_cast<size_t>(hardScattering.NVariables()));
    achilles::Vegas vegas(map, node["Vegas"]);
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
