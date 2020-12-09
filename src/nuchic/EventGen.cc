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
    cascade = std::make_unique<Cascade>(config["Cascade"].as<Cascade>());

    auto scatteringNode = config["Main"]["Hard Scattering"];
    auto runMode = config["Main"]["Run Mode"].as<nuchic::RunMode>();
    scattering = HardScatteringFactory::Create(scatteringNode["Model"].as<std::string>(),
            scatteringNode, runMode);
    if(runMode == RunMode::FixedAngle)
        scattering -> SetScatteringAngle(config["Main"]["Angle"].as<double>());

    // Setup Vegas
    nuchic::AdaptiveMap map(static_cast<size_t>(scattering->NVariables() + beam->NVariables()));
    integrator = Vegas(map, config["Initialize"]);

    auto output = config["Main"]["Output"];
    if(output["Format"].as<std::string>() == "Nuchic") {
        writer = std::make_unique<NuchicWriter>(output["Name"].as<std::string>());
    }
    writer -> WriteHeader(configFile);

    hist = Histogram(2, -0.5, 1.5, "xsec");
    rng = std::make_shared<randutils::mt19937_rng>();
}

void nuchic::EventGen::Initialize() {
    auto func = [&](const std::vector<double> &x, const double &wgt) {
        return Calculate(x, wgt);
    };
    integrator(func);
}

void nuchic::EventGen::GenerateEvents() {
    integrator.Clear();
    integrator.Set(config["EventGen"]);
    outputEvents = true;
    runCascade = config["Cascade"]["Run"].as<bool>();

    auto func = [&](const std::vector<double> &x, const double &wgt) {
        auto niterations = config["EventGen"]["iterations"].as<double>();
        return Calculate(x, wgt/niterations);
    };
    integrator(func);

    hist.Save("multi");
}

double nuchic::EventGen::Calculate(const std::vector<double> &rans, const double &wgt) {
    // Initialize the event, which generates the nuclear configuration
    // and initializes the beam particle for the event
    std::vector<double> beamRans(rans.begin(), rans.begin() + beam -> NVariables());
    Event event(nucleus, beam, beamRans, wgt);

    // Generate phase space
    scattering -> GeneratePhaseSpace(rans, event);
    if(event.PhaseSpace().weight == 0) return 0;

    // Calculate the hard cross sections and select one for initial state
    scattering -> CrossSection(event);
    if(!scattering -> InitializeEvent(event))
        return 0;

    // Run the cascade if needed
    if(runCascade) {
        cascade -> Evolve(event);
    } else {
        for(auto & nucleon : event.CurrentNucleus()->Nucleons()) {
            if(nucleon.Status() == ParticleStatus::propagating) {
                nucleon.SetStatus(ParticleStatus::escaped);
            }
        }
    }

    // Write out events
    if(outputEvents) {
        writer -> Write(event);
    }

    return event.Weight()/wgt;
}

// Move this to hardscattering starting here:
//    double totalXSec = std::accumulate(xsecwgt.begin(), xsecwgt.end(), 0.0);
//
//    for particle in particles:
//        if momentum[2] (initial state nucleon) < nucleus.kf(particle.position):
//            if proton:
//                add proton to list of allowed proton initial states
//            if neutron:
//                add neutron to list of allowed neutron initial states
//
//    totalXSec = (n_protons_allowed*rho*xsecwgt[0] + n_neutrons_allowed*xsecwgt[1]);
//    // two body
//    totalXSec = (n_protons_allowed*(n_protons_allowed-1)*xsecwgt[0]
//                 + n_protons_allowed*n_neutrons_allowed*xsecwgt[1]
//                 + n_neutrons_allowed*(n_neutrons_allowed-1)*xsecwgt[2]);
//    // Select struck nucleon(s)
//    bool proton = rng -> uniform(0.0, totalXSec) < xsecwgt[0];
//    if(proton) {
//        select from list of allowed initial states
//
//        auto idx = rng -> pick(nucleus -> ProtonsIDs());
//        particles[idx].SetStatus(ParticleStatus::initial_state);
//        particles[idx].SetMomentum(psPoint.momentum[2]);
//        Particle outNucl(particles[idx]);
//        outNucl.SetStatus(ParticleStatus::propagating);
//        outNucl.SetMomentum(psPoint.momentum[3]);
//        particles.push_back(outNucl);
//    } else {
//        auto idx = rng -> pick(nucleus -> NeutronsIDs());
//        particles[idx].SetStatus(ParticleStatus::initial_state);
//        particles[idx].SetMomentum(psPoint.momentum[2]);
//        Particle outNucl(particles[idx]);
//        outNucl.SetStatus(ParticleStatus::propagating);
//        outNucl.SetMomentum(psPoint.momentum[3]);
//        particles.push_back(outNucl);
//    }
// To here
