#include "nuchic/EventGen.hh"
#include "nuchic/Event.hh"
#include "nuchic/EventWriter.hh"
#include "nuchic/HardScatteringFactory.hh"
#include "nuchic/HardScattering.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Beams.hh"
#include "nuchic/Cascade.hh"
#include "nuchic/Particle.hh"

#include "yaml-cpp/yaml.h"

nuchic::EventGen::EventGen(const std::string &configFile) : runCascade{false}, outputEvents{false} {
    config = YAML::LoadFile(configFile);

    beam = std::make_shared<Beam>(config["Beams"].as<Beam>());
    nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());
    // cascade = std::make_unique<Cascade>(config["Cascade"].as<Cascade>());

    auto scatteringNode = config["Main"]["Hard Scattering"];
    auto runMode = config["Main"]["Run Mode"].as<nuchic::RunMode>();
    scattering = HardScatteringFactory::Create(scatteringNode["Model"].as<std::string>(),
            scatteringNode, beam, nucleus, runMode);
    if(runMode == RunMode::FixedAngle)
        scattering -> SetScatteringAngle(config["Main"]["Angle"].as<double>());

    // Setup Vegas
    nuchic::AdaptiveMap map(static_cast<size_t>(scattering->NVariables()));
    integrator = Vegas(map, config["Initialize"]);

    auto output = config["Main"]["Output"];
    if(output["Format"].as<std::string>() == "Nuchic") {
        writer = std::make_unique<NuchicWriter>(output["Name"].as<std::string>());
    }
    writer -> WriteHeader(configFile);
}

void nuchic::EventGen::Initialize() {
    auto func = [&](const std::vector<double> &x, const double wgt) {
        return Calculate(x, wgt);
    };
    integrator(func);
}

void nuchic::EventGen::GenerateEvents() {
    integrator.Clear();
    integrator.Set(config["EventGen"]);
    outputEvents = true;
    runCascade = config["Cascade"]["Run"].as<bool>();

    auto func = [&](const std::vector<double> &x, const double wgt) {
        auto niterations = config["EventGen"]["iterations"].as<double>();
        return Calculate(x, wgt/niterations);
    };
    integrator(func);
}

double nuchic::EventGen::Calculate(const std::vector<double> &x, const double &/*wgt*/) {
    // TODO: Choose nucleon here!!!
    static constexpr double conv = 1e6;
    auto particles = scattering -> GeneratePhaseSpace(x);
    if(particles[2].E() <= 0) return 0.0;
    double pswgt = scattering -> PhaseSpaceWeight(particles);
    if(pswgt == 0) return pswgt;
    double xsecwgt = scattering -> CrossSection(particles);

    // TODO: Figure out how to properly give the right nucleon a kick
    if(runCascade) {
        // Select a random nucleon to kick 
    }
    // TODO: write out events to file
    if(outputEvents) {
        Event event(std::move(particles), pswgt*xsecwgt*conv);
        writer -> Write(event);
    }

    return pswgt*xsecwgt*conv;
}
