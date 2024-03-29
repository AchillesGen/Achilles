#include "Achilles/EventGen.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/HardScatteringFactory.hh"
#include "Achilles/HardScattering.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Units.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/ComplexFmt.hh"
#include "Achilles/Units.hh"
#include "Achilles/Channels.hh"

#ifdef ENABLE_BSM
#include "plugins/Sherpa/SherpaInterface.hh"
#endif

#ifdef ENABLE_HEPMC3
#include "plugins/HepMC3/HepMC3EventWriter.hh"
#include "plugins/NuHepMC/NuHepMCWriter.hh"
#endif

#include "yaml-cpp/yaml.h"


achilles::EventGen::EventGen(const std::string &configFile,
                             std::vector<std::string> shargs) : config{configFile} {

    // Turning off decays in Sherpa. This is a temporary fix until we can get ISR and FSR properly working in SHERPA.
    runDecays = false;

    // Setup random number generator
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config.Exists("Options/Initialize/Seed"))
        if(config.GetAs<int>("Options/Initialize/Seed") > 0)
            seed = config.GetAs<unsigned int>("Options/Initialize/Seed");
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Setup unweighter
    unweighter = UnweighterFactory::Initialize(config.GetAs<std::string>("Options/Unweighting/Name"),
                                               config["Options/Unweighting"]);

    // Load initial state, massess
    spdlog::trace("Initializing the beams");
    beam = std::make_shared<Beam>(config.GetAs<Beam>("Beams"));
    nucleus = std::make_shared<Nucleus>(config.GetAs<Nucleus>("Nucleus"));

    // Set potential for the nucleus
    auto potential_name = config.GetAs<std::string>("Nucleus/Potential/Name");
    auto potential = achilles::PotentialFactory::Initialize(potential_name,
                                                          nucleus,
                                                          config["Nucleus/Potential"]);
    nucleus -> SetPotential(std::move(potential));
    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config.GetAs<bool>("Cascade/Run"));
    if(config.GetAs<bool>("Cascade/Run")) {
        cascade = std::make_unique<Cascade>(config.GetAs<Cascade>("Cascade"));
    } else {
        cascade = nullptr;
    }

    // Initialize the lepton final states
    spdlog::debug("Initializing the leptonic final states");
    auto leptonicProcess = config.GetAs<achilles::Process_Info>("Process");
    // TODO: Handle the beam initial state better
    if(beam -> BeamIDs().size() > 1)
        throw std::runtime_error("Multiple processes are not implemented yet. Please use only one beam.");
    leptonicProcess.m_ids.insert(leptonicProcess.m_ids.begin(),
                                 beam -> BeamIDs().begin(), beam -> BeamIDs().end());

    // Initialize the nuclear model
    spdlog::debug("Initializing nuclear model");
    const auto model_name = config.GetAs<std::string>("NuclearModel/Model");
    auto nuclear_model = NuclearModelFactory::Initialize(model_name, config.Root());
    nuclear_model -> AllowedStates(leptonicProcess);
    spdlog::debug("Process: {}", leptonicProcess);

#ifdef ENABLE_BSM
    // Initialize sherpa processes
    p_sherpa = new achilles::SherpaInterface();
    std::string model = config.GetAs<std::string>("Process/Model");
    std::string param_card = config.GetAs<std::string>("Process/ParamCard");
    int qed = 0;
    if(config.Exists("Process/QEDShower"))
        if(config.GetAs<bool>("Process/QEDShower"))
            qed = 3;
    shargs.push_back(fmt::format("CSS_EW_MODE={}", qed));
    if(model == "SM") model = "SM_Nuc";
    shargs.push_back("MODEL=" + model);
    shargs.push_back("UFO_PARAM_CARD=" + param_card);
    shargs.push_back(fmt::format("BEAM_2={}", 11));
    shargs.push_back(fmt::format("BEAM_ENERGY_2={}", 20)); 
    p_sherpa -> Initialize(shargs);
    spdlog::debug("Initializing leptonic currents");
    if(!p_sherpa -> InitializeProcess(leptonicProcess)) {
        spdlog::error("Cannot initialize hard process");
        exit(1);
    }
    leptonicProcess.m_mom_map = p_sherpa -> MomentumMap(leptonicProcess.Ids());
#else
    // Dummy call to remove unused error
    (void)shargs.size();
    leptonicProcess.m_mom_map[0] = leptonicProcess.Ids()[0];
    leptonicProcess.m_mom_map[1] = leptonicProcess.Ids()[1];
    leptonicProcess.m_mom_map[2] = leptonicProcess.Ids()[2];
    leptonicProcess.m_mom_map[3] = leptonicProcess.Ids()[3];
#endif

    // Initialize hard cross-sections
    spdlog::debug("Initializing hard interaction");
    scattering = std::make_shared<HardScattering>();
    scattering -> SetProcess(leptonicProcess);
#ifdef ENABLE_BSM
    scattering -> SetSherpa(p_sherpa);
#endif
    scattering -> SetNuclear(std::move(nuclear_model));

    // Setup channels
    spdlog::debug("Initializing phase space");
    std::vector<double> masses = scattering -> Process().Masses();
    spdlog::trace("Masses = [{}]", fmt::join(masses.begin(), masses.end(), ", "));
    if(config.Exists("TestingPS")) {
        Channel<FourVector> channel = BuildChannelTest(config["TestingPS"], beam);
        integrand.AddChannel(std::move(channel));
    } else {
#ifndef ENABLE_BSM
        if(scattering -> Process().Multiplicity() == 4) {
            Channel<FourVector> channel0 = BuildChannel<TwoBodyMapper>(scattering -> Nuclear(), 2, 2,
                                                                       beam, masses);
            integrand.AddChannel(std::move(channel0));
        } else {
            const std::string error = fmt::format("Leptonic Tensor can only handle 2->2 processes without "
                                                  "BSM being enabled. "
                                                  "Got a 2->{} process", leptonicProcess.m_ids.size());
            throw std::runtime_error(error);
        }
#else
        auto channels = p_sherpa -> GenerateChannels(scattering -> Process().Ids());
        size_t count = 0;
        for(auto & chan : channels) {
            Channel<FourVector> channel = BuildGenChannel(scattering -> Nuclear(), 
                                                          scattering -> Process().m_ids.size(), 2,
                                                          beam, std::move(chan), masses);
            integrand.AddChannel(std::move(channel));
            spdlog::info("Adding Channel{}", count++);
        }
#endif
    }

    // Setup Multichannel integrator
    // auto params = config["Integration"]["Params"].as<MultiChannelParams>();
    integrator = MultiChannel(integrand.NDims(), integrand.NChannels(), {1000, 2});

    // Decide whether to rotate events to be measured w.r.t. the lepton plane
    if(config.Exists("Main/DoRotate"))
        doRotate = config.GetAs<bool>("Main/DoRotate");

    // Setup Cuts
    if(config.Exists("Main/HardCuts")) {
        doHardCuts = config.GetAs<bool>("Main/HardCuts");
        spdlog::info("Apply hard cuts? {}", doHardCuts);
        hard_cuts = config.GetAs<achilles::CutCollection>("HardCuts");
    }

    // if(config["Main"]["EventCuts"])
    //     doEventCuts = config["Main"]["EventCuts"].as<bool>();
    // spdlog::info("Apply event cuts? {}", doEventCuts);
    // event_cuts = config["EventCuts"].as<achilles::CutCollection>();

    // Setup outputs
    bool zipped = true;
    if(config.Exists("Main/Output/Zipped"))
        zipped = config.GetAs<bool>("Main/Output/Zipped");
    auto format = config.GetAs<std::string>("Main/Output/Format");
    auto name = config.GetAs<std::string>("Main/Output/Name");
    spdlog::trace("Outputing as {} format", format);
    if(format == "Achilles") {
        writer = std::make_unique<AchillesWriter>(name, zipped);
#ifdef ENABLE_HEPMC3
    } else if(format == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(name, zipped);
    } else if(format == "NuHepMC") {
        writer = std::make_unique<NuHepMCWriter>(name, zipped);
#endif
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}", format);
        throw std::runtime_error(msg);
    }
    writer -> WriteHeader(configFile);
}

void achilles::EventGen::Initialize() {
    // TODO: Clean up loading of previous results
    auto func = [&](const std::vector<FourVector> &mom, const double &wgt) {
        return GenerateEvent(mom, wgt);
    };
    // TODO: Loading the saved data is broken
    // try {
    //     YAML::Node old_results = YAML::LoadFile("results.yml");
    //     integrator = old_results["Multichannel"].as<MultiChannel>();
    //     integrand = old_results["Channels"].as<Integrand<FourVector>>();
    //     YAML::Node results;
    //     results["Multichannel"] = integrator;
    //     results["Channels"] = integrand;
    //     integrand.Function() = func;
    // } catch(const YAML::BadFile &e) {
        spdlog::info("Initializing integrator.");
        integrand.Function() = func;
        if(config.Exists("Initialize/Accuracy"))
            integrator.Parameters().rtol = config["Initialize"]["Accuracy"].as<double>();
        integrator.Optimize(integrand);
        integrator.Summary();

        YAML::Node results;
        results["Multichannel"] = integrator;
        results["Channels"] = integrand;

        std::ofstream fresults("results.yml");
        fresults << results;
        fresults.close();
    // }
}

void achilles::EventGen::GenerateEvents() {
    outputEvents = true;
    runCascade = config["Cascade/Run"].as<bool>();
    integrator.Parameters().ncalls = config["Main/NEvents"].as<size_t>();
    integrator(integrand);
    fmt::print("\n");
    auto result = integrator.Summary();
    fmt::print("Integral = {:^8.5e} +/- {:^8.5e} ({:^8.5e} %)\n",
               result.results.back().Mean(), result.results.back().Error(),
               result.results.back().Error() / result.results.back().Mean()*100);
    fmt::print("Unweighting efficiency: {:^8.5e} %\n",
               unweighter->Efficiency() * 100);
}

double achilles::EventGen::GenerateEvent(const std::vector<FourVector> &mom, const double &wgt) {
    if(outputEvents) {
        static constexpr size_t statusUpdate = 1000;
        if(unweighter->Accepted() % statusUpdate == 0) {
            fmt::print("Generated {} / {} events\r",
                       unweighter->Accepted(),
                       config["Main/NEvents"].as<size_t>());
        }
    }
    // Initialize the event, which generates the nuclear configuration
    // and initializes the beam particle for the event
    Event event(nucleus, mom, wgt);

    // Initialize the particle ids for the processes
    const auto pids = scattering -> Process().m_ids;

    // Setup flux value
    event.Flux() = beam -> EvaluateFlux(pids[0], mom[1]);

    spdlog::trace("Event Phase Space:");
    size_t idx = 0;
    for(const auto &momentum : event.Momentum()) {
        spdlog::trace("\t{}: {} (M2 = {})", ++idx, momentum, momentum.M2());
    }

    // Calculate the hard cross sections and select one for initial state
    spdlog::trace("Calculating cross section");

    // Obtain the cross section for the event 
    auto xsecs = scattering -> CrossSection(event);

    // Initialize the event
    spdlog::trace("Filling the event");
    if(!scattering -> FillEvent(event, xsecs)) {
        if(outputEvents) {
            event.SetMEWeight(0);
            event.CalcWeight();
            spdlog::trace("Outputting the event");
            writer -> Write(event);

            // Update number of calls needed to ensure the number of generated events
            // is the same as that requested by the user
            integrator.Parameters().ncalls++;
        }

#ifdef ENABLE_BSM
        p_sherpa->Reset();
#endif
        return 0;
    }

#ifdef ACHILLES_EVENT_DETAILS
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
#endif

    event.CalcWeight();
    spdlog::trace("Weight: {}", event.Weight());

    // Perform hard cuts
    if(doHardCuts) {
        spdlog::trace("Making hard cuts");
        if(!MakeCuts(event)) {
            // Short-circuit the evaluation
            // We want Vegas to adapt to avoid these points, i.e.,
            // the integrand should be interpreted as zero in this region
            if(outputEvents) {
                event.SetMEWeight(0);
                event.CalcWeight();
                writer -> Write(event);
                // Update number of calls needed to ensure the number of generated events
                // is the same as that requested by the user
                integrator.Parameters().ncalls++;
            }

#ifdef ENABLE_BSM
            p_sherpa->Reset();
#endif
            return 0;
        }
    }

    // TODO: Move to after unweighting?
    // Run the cascade if needed
    if(runCascade) {
#ifdef ACHILLES_EVENT_DETAILS
        spdlog::trace("Hadrons:");
        idx = 0;
        for(const auto &particle : event.Hadrons()) {
            if(particle.Status() == ParticleStatus::initial_state
                && particle.ID() == PID::proton())
                spdlog::trace("\t{}: {}", idx, particle);
            ++idx;
        }
#endif
        spdlog::trace("Runnning cascade");
        cascade -> Evolve(&event);

#ifdef ACHILLES_EVENT_DETAILS
        spdlog::trace("Hadrons (Post Cascade):");
        idx = 0;
        for(const auto &particle : event.Hadrons()) {
            spdlog::trace("\t{}: {}", ++idx, particle);
        }
#endif
    } else {
        for(auto & nucleon : event.CurrentNucleus()->Nucleons()) {
            if(nucleon.Status() == ParticleStatus::propagating) {
                nucleon.Status() = ParticleStatus::final_state;
            }
        }
    }

    // Write out events
    if(outputEvents) {
        // TODO: Handle MEC case
        // Setup target nucleus in history
        auto init_nuc = event.CurrentNucleus()->InitParticle();
        Particle init_had;
        for(const auto &nucleon : event.CurrentNucleus()->Nucleons()) {
            if(nucleon.Status() == ParticleStatus::initial_state) {
                init_had = nucleon;
                break;
            }
        }
        event.History().AddVertex(init_had.Position(), {init_nuc}, {init_had}, EventHistory::StatusCode::target);
        // Setup beam in history
        auto init_lep = event.Leptons()[0];
        auto init_beam = init_lep;
        init_beam.Status() = ParticleStatus::beam;
        const double max_energy = beam->MaxEnergy();
        init_beam.Momentum() = {max_energy, 0, 0, max_energy};
        event.History().AddVertex({}, {init_beam}, {init_lep}, EventHistory::StatusCode::beam);
#ifdef ENABLE_BSM
        // Running Sherpa interface if requested
        // Only needed when generating events and not optimizing the multichannel
        // TODO: Move to after unweighting?
        if(event.Weight() > 0) {
            if(runDecays)
                p_sherpa -> GenerateEvent(event);
            else {
                // TODO: Properly build history including the cascade
                std::vector<Particle> final;
                for(const auto &part : event.Particles()) {
                    if(part.IsFinal()) final.push_back(part);
                }
                event.History().AddVertex(init_had.Position(), {init_had, init_lep}, {final},
                                          EventHistory::StatusCode::primary); 
            }
        }
#else
        if(event.Weight() > 0) {
            // TODO: Properly build history including the cascade
            std::vector<Particle> final;
            for(const auto &part : event.Particles()) {
                if(part.IsFinal()) final.push_back(part);
            }
            event.History().AddVertex(init_had.Position(), {init_had, init_lep}, {final},
                                      EventHistory::StatusCode::primary); 
        }
#endif
        // TODO: Get remnant working
        // Setup remnant in history
        // auto recoilMom = init_nuc.Momentum();
        // for(size_t i = 1; i < event.Leptons().size(); ++i) {
        //     recoilMom -= event.Leptons()[i].Momentum();
        // }
        // for(size_t i = 1; i < event.Hadrons().size(); ++i) {
        //     recoilMom -= event.Hadrons()[i].Momentum();
        // }
        // auto remnant = Particle(event.Remnant().PID(), recoilMom); 
        // event.History().Primary()->AddOutgoing(remnant);

        // Rotate cuts into plane of outgoing electron before writing
        if (doRotate)
            Rotate(event);
        // Perform event-level final cuts before writing
        bool outputCurrentEvent = true;
        // if(doEventCuts){
        //     spdlog::debug("Making event cuts");
        //     outputCurrentEvent = MakeEventCuts(event);
        // }

        if(outputCurrentEvent) {
            // Keep a running total of the number of surviving events
            event.Finalize();

            if(!unweighter->AcceptEvent(event)) {
                // Update number of calls needed to ensure the number of generated events
                // is the same as that requested by the user
                integrator.Parameters().ncalls++;
            }
            writer -> Write(event);
        }
    } else {
        unweighter->AddEvent(event);
    }

#ifdef ENABLE_BSM
    p_sherpa->Reset();
#endif

    // Always return the weight when the event passes the initial hard cut.
    // Even if events do not survive the final event-level cuts, Vegas should
    // still interpret the integrand as nonzero in this region.
    return event.Weight();
}

bool achilles::EventGen::MakeCuts(Event &event) {
    return hard_cuts.EvaluateCuts(event.Particles());
}

// TODO: Create Analysis level cuts
/*
bool achilles::EventGen::MakeEventCuts(Event &event) {
    // Run through all particles in the event
    for (const auto& pair : event_cuts) {
        auto pid = pair.first;
        auto cut = pair.second;
        bool pid_passed = false;
        for (const auto& particle : event.Particles()){
            // Restrict to matching final-state particles
            if(particle.IsFinal() && particle.ID() == pid)
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
}*/

void achilles::EventGen::Rotate(Event &event) {
    // Isolate the azimuthal angle of the outgoing electron
    double phi = 0.0;
    for(const auto & particle : event.Particles()){
        if(particle.ID() == PID::electron() && particle.IsFinal()){
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
