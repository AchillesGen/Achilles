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
    // Setup random number generator
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["seed"])
        seed = config["Initialize"]["seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load initial states
    beam = std::make_shared<Beam>(config["Beams"].as<Beam>());
    nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());

    // Set potential for the nucleus
    auto potential_name = config["Nucleus"]["Potential"]["Name"].as<std::string>();
    auto potential = nuchic::PotentialFactory::Initialize(potential_name,
                                                          nucleus,
                                                          config["Nucleus"]["Potential"]);
    nucleus -> SetPotential(std::move(potential));

    // Event counts
    total_events = config["EventGen"]["TotalEvents"].as<size_t>();
    nevents = 0; // Initialize to zero
    max_batch = config["EventGen"]["MaxBatch"].as<size_t>();

    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Run"].as<bool>());
    if(config["Cascade"]["Run"].as<bool>()) {
        cascade = std::make_unique<Cascade>(config["Cascade"].as<Cascade>());
    } else {
        cascade = nullptr;
    }

    // Initialize hard cross-sections
    auto scatteringNode = config["Main"]["Hard Scattering"];
    auto runMode = config["Main"]["Run Mode"].as<nuchic::RunMode>();
    scattering = HardScatteringFactory::Create(scatteringNode["Model"].as<std::string>(),
            scatteringNode, runMode);
    if(runMode == RunMode::FixedAngle)
        scattering -> SetScatteringAngle(config["Main"]["Angle"].as<double>()*1.0_deg);
    if(runMode == RunMode::FixedAngleEnergy) {
        scattering -> SetScatteringAngle(config["Main"]["Angle"].as<double>()*1.0_deg);
        scattering -> SetFinalLeptonEnergy(config["Main"]["ELepFinal"].as<double>());
    }

    // Setup Vegas
    nuchic::AdaptiveMap map(static_cast<size_t>(scattering->NVariables() + beam->NVariables()));
    integrator = Vegas(map, config["Initialize"]);

    // Decide whether to rotate events to be measured w.r.t. the lepton plane
    if(config["Main"]["DoRotate"])
        doRotate = config["Main"]["DoRotate"].as<bool>();

    // Setup Cuts
    if(config["Main"]["HardCuts"])
        doHardCuts = config["Main"]["HardCuts"].as<bool>();
    spdlog::info("Apply hard cuts? {}", doHardCuts);
    hard_cuts = config["HardCuts"].as<nuchic::Cuts>();

    if(config["Main"]["EventCuts"])
        doEventCuts = config["Main"]["EventCuts"].as<bool>();
    spdlog::info("Apply event cuts? {}", doEventCuts);
    event_cuts = config["EventCuts"].as<nuchic::Cuts>();

    // Setup outputs
    auto output = config["Main"]["Output"];
    if(output["Format"].as<std::string>() == "Nuchic") {
        bool zipped = true;
        if(output["Zipped"])
            zipped = output["Zipped"].as<bool>();
        writer = std::make_unique<NuchicWriter>(output["Name"].as<std::string>(), zipped);
    }
    writer -> WriteHeader(configFile);

    hist = Histogram(1000, 0.0, 1000.0, "xsec");
    hist2 = Histogram(100, 0.0, 1400.0, "tpe");
    hist3 = Histogram(13, -0.5, 12.5, "np");
    hist4 = Histogram(100, 0.0, 1400.0, "erec");
    
}

void nuchic::EventGen::Initialize() {
    spdlog::info("Initializing vegas integrator.");
    size_t batch_count = 0;
    auto func = [&](const std::vector<double> &x, const double &wgt) {
        return Calculate(x, wgt, batch_count);
    };
    integrator(func);
}

void nuchic::EventGen::GenerateEvents() {
    // integrator.Clear();
    integrator.Set(config["EventGen"]);
    outputEvents = true;
    runCascade = config["Cascade"]["Run"].as<bool>();
    spdlog::info("Starting generating of n >= {} total events", total_events);
    spdlog::info("Using a maximum of {} total Vegas batches.", max_batch);
    // Run integrator in batches until the desired number of events are found
    size_t batch_count = 1;
    while ((nevents < total_events) & (batch_count <= max_batch)){
        integrator.Clear();  // Reset integrator for each batch
        auto func = [&](const std::vector<double> &x, const double &wgt) {
            auto niterations = config["EventGen"]["iterations"].as<double>();
            return Calculate(x, wgt/niterations, batch_count);
        };
        spdlog::info("Running vegas batch number {}", batch_count);
        integrator(func);
        spdlog::info("Total events so far: {}/{}", nevents, total_events);
        batch_count += 1;
    }
    if (batch_count >= max_batch){
        spdlog::info("Stopping after reaching max batch threshold.");
    }

    hist.Save("multi");
    hist2.Save("tpe");
    hist3.Save("np");
    hist4.Save("erec");
    
}

double nuchic::EventGen::Calculate(const std::vector<double> &rans, const double &wgt, const size_t &batch) {
    // Initialize the event, which generates the nuclear configuration
    // and initializes the beam particle for the event
    std::vector<double> beamRans(rans.begin(), rans.begin() + beam -> NVariables());
    Event event(nucleus, beam, beamRans, wgt);
    event.SetBatch(batch);

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

    // Perform hard cuts
    if(doHardCuts) {
        spdlog::debug("Making hard cuts");
        if(!MakeCuts(event))
            // Short-circuit the evaluation
            // We want Vegas to adapt to avoid these points, i.e.,
            // the integrand should be interpreted as zero in this region
            return 0;
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

    // Write out events
    if(outputEvents) {
        // Rotate cuts into plane of outgoing electron before writing
        if (doRotate)
            Rotate(event);

        bool writeEvent = true;
        // Perform event-level final cuts before writing
        if(doEventCuts){
            spdlog::debug("Making event cuts");
            writeEvent = MakeEventCuts(event);
        }

        if(writeEvent) {
            // Keep a running total of the number of surviving events
            nevents += 1;
            spdlog::debug("Found event: {}/{}", nevents, total_events);
            event.Finalize();
            writer -> Write(event);
            const auto omega = event.Leptons()[0].E() - event.Leptons()[1].E();
            hist.Fill(omega, event.Weight()/(2*M_PI));
	    const auto eps=21.0;
	    const auto erec= (Constant::mN*eps + Constant::mN*event.Leptons()[1].E())
		    /(Constant::mN-event.Leptons()[1].E()+event.Leptons()[1].Pz() );
	    const auto cosAngle=event.Leptons()[1].Pz()/event.Leptons()[1].E();
	    auto theta = std::acos(cosAngle)*180/M_PI;
	    if (theta > 17+ 7_GeV/event.Leptons()[1].E() && event.Leptons()[1].E()>400) { 
            	hist4.Fill(erec, event.Weight());
	    }
	    
	    int nprotons = 0;
	    double tpe= 0;
	    for(const auto &nucleon : event.CurrentNucleus()->Nucleons()) {
		    if(nucleon.Status() == ParticleStatus::escaped && nucleon.ID() == PID::proton()){
                            double kin = nucleon.Momentum().E() - Constant::mN+eps +event.Leptons()[1].E();
			    if (kin > tpe) tpe=kin;
			    nprotons++;
		    }
	    }
            if (nprotons==1) {	    	    
	          hist2.Fill(tpe, event.Weight());
	    }	  
	    hist3.Fill(nprotons, event.Weight());
	    
        }
    }

    // Always return the weight when the event passes the initial hard cut.
    // Even if events do not survive the final event-level cuts, Vegas should
    // still interpret the integrand as nonzero in this region.
    return event.Weight()/wgt;
}

bool nuchic::EventGen::MakeCuts(Event &event) {
    // Run through all particles in the event
    for(const auto &particle : event.Particles())
        // Only apply cuts to final-state particles
        if(particle.IsFinal())
            if(hard_cuts.find(particle.ID()) != hard_cuts.end()){
                // Reject the event if a single particle fails a cut
                if(!hard_cuts[particle.ID()](particle.Momentum())){
                    return false;
                }
            }
    return true;
}

bool nuchic::EventGen::MakeEventCuts(Event &event) {
    // Run through all particles in the event
    for (const auto& pair : event_cuts) {
        auto pid = pair.first;
        auto cut = pair.second;
        bool pid_passed = false;
        for (const auto& particle : event.Particles()){
            // Restrict to matching final-state particles
            if(particle.IsFinal() & (particle.ID() == pid))
                // Keep: at least one particle (of a given PID) survives the cut
                if(cut(particle.Momentum())){
                    pid_passed = true;
                    break;
                }
        }
        // Reject: no particles (of a given PID) satisfy the cut
        if(!pid_passed)
            return false;
    }
    return true;
}




void nuchic::EventGen::Rotate(Event &event) {
    // Isolate the azimuthal angle of the outgoing electron
    double phi = 0.0;
    for(const auto & particle : event.Particles()){
        if((int(particle.ID()) == 11) & particle.IsFinal()){
            phi = particle.Momentum().Phi();
        }
    }
    // Rotate the coordiantes of particles so that all azimuthal angles phi are
    // measured with respect to the leptonic plane
    std::array<double, 9> rotation = {
        cos(phi),  sin(phi), 0,
        -sin(phi), cos(phi), 0,
        0,         0,        1};
    event.Rotate(rotation);
}
