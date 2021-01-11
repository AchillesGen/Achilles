#include "nuchic/EventGen.hh"
#include "nuchic/Event.hh"
#include "nuchic/EventWriter.hh"
#include "nuchic/HardScatteringFactory.hh"
#include "nuchic/HardScattering.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Beams.hh"
#include "nuchic/Cascade.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Units.hh"

#include "yaml-cpp/yaml.h"

nuchic::EventGen::EventGen(const std::string &configFile) : runCascade{false}, outputEvents{false} {
    config = YAML::LoadFile(configFile);

    // Load initial states
    beam = std::make_shared<Beam>(config["Beams"].as<Beam>());
    nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());

    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Run"].as<bool>());
    if(config["Cascade"]["Run"].as<bool>()) {
        cascade = std::make_unique<Cascade>(config["Cascade"].as<Cascade>());
    } else {
        cascade = nullptr;
    }

    // Setup random number generator
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["Seed"])
        seed = config["Initialize"]["Seed"].as<unsigned int>();
    rng = std::make_shared<randutils::mt19937_rng>(seed);

    // Initialize hard cross-sections
    auto scatteringNode = config["Main"]["Hard Scattering"];
    auto runMode = config["Main"]["Run Mode"].as<nuchic::RunMode>();
    scattering = HardScatteringFactory::Create(scatteringNode["Model"].as<std::string>(),
            scatteringNode, runMode, rng);
    if(runMode == RunMode::FixedAngle)
        scattering -> SetScatteringAngle(config["Main"]["Angle"].as<double>()*1.0_deg);

    // Setup Vegas
    nuchic::AdaptiveMap map(static_cast<size_t>(scattering->NVariables() + beam->NVariables()));
    integrator = Vegas(map, config["Initialize"], rng);

    // Setup Cuts
    doCuts = config["Main"]["DoCuts"].as<bool>();
    cuts = config["Cuts"].as<nuchic::Cuts>();

    // Setup outputs
    auto output = config["Main"]["Output"];
    if(output["Format"].as<std::string>() == "Nuchic") {
        writer = std::make_unique<NuchicWriter>(output["Name"].as<std::string>());
    }
    writer -> WriteHeader(configFile);

    hist = Histogram(1000, 0.0, 1000.0, "xsec");
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
    spdlog::debug("Generating phase space");
    scattering -> GeneratePhaseSpace(rans, event);
    if(event.PhaseSpace().weight == 0) return 0;

    // Calculate the hard cross sections and select one for initial state
    spdlog::debug("Calculating cross section");
    scattering -> CrossSection(event);
    if(!scattering -> InitializeEvent(event))
        return 0;

    spdlog::trace("Event Phase Space:");
    size_t idx = 0;
    for(const auto &mom : event.PhaseSpace().momentum) {
        spdlog::trace("\t{}: {}", ++idx, mom);
    }

    spdlog::trace("Leptons:");
    idx = 0;
    for(const auto &particle : event.Leptons()) {
        spdlog::trace("\t{}: {}", ++idx, particle);
    }

    spdlog::trace("Hadrons:");
    idx = 0;
    for(const auto &particle : event.Hadrons()) {
        spdlog::trace("\t{}: {}", ++idx, particle);
    }

    // Run the cascade if needed
    if(runCascade) {
        spdlog::debug("Runnning cascade");
        cascade -> Evolve(&event);
    } else {
        for(auto & nucleon : event.CurrentNucleus()->Nucleons()) {
            if(nucleon.Status() == ParticleStatus::propagating) {
                nucleon.Status() = ParticleStatus::escaped;
            }
        }
    }

    // Preform cuts
    if(doCuts) {
        spdlog::debug("Making cuts");
        if(!MakeCuts(event))
            return 0;
    }

    // Write out events
    if(outputEvents) {
        event.Finalize();
        writer -> Write(event);
        const auto omega = event.Leptons()[0].E() - event.Leptons()[1].E();
        hist.Fill(omega, event.Weight()/(2*M_PI));
	
    }

    return event.Weight()/wgt;
}

bool nuchic::EventGen::MakeCuts(Event &event) {
    for(const auto &particle : event.Particles())
        if(particle.IsFinal())
            if(cuts.find(particle.ID()) != cuts.end())
                if(!cuts[particle.ID()](particle.Momentum()))
                    return false;
    return true;
}
