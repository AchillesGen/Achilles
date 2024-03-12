#include "Achilles/EventGen.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/HardScatteringFactory.hh"
#include "Achilles/HardScattering.hh"
#include "Achilles/Logging.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Units.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/ComplexFmt.hh"
#include "Achilles/Units.hh"

// TODO: Turn this into a factory to reduce the number of includes
#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/BeamMapper.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/PhaseSpaceMapper.hh"
#include "Achilles/QuasielasticTestMapper.hh"

#ifdef ENABLE_BSM
#include "plugins/Sherpa/Channels1.hh"
#include "plugins/Sherpa/Channels3.hh"
#include "plugins/Sherpa/SherpaMEs.hh"
#endif

#ifdef ENABLE_HEPMC3
#include "plugins/HepMC3/HepMC3EventWriter.hh"
#endif

#include "yaml-cpp/yaml.h"

achilles::Channel<achilles::FourVector> BuildChannelTest(const YAML::Node &node, std::shared_ptr<achilles::Beam> beam) {
    achilles::Channel<achilles::FourVector> channel;
    channel.mapping = std::make_unique<achilles::QuasielasticTestMapper>(node, beam);
    achilles::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = achilles::Vegas(map, {});
    return channel;
}

template<typename T>
achilles::Channel<achilles::FourVector> BuildChannel(achilles::NuclearModel *model, size_t nlep, size_t nhad,
                                                     std::shared_ptr<achilles::Beam> beam,
                                                     const std::vector<double> &masses,
                                                     std::optional<double> gauge_boson_mass = std::nullopt) {
    achilles::Channel<achilles::FourVector> channel;
    channel.mapping = achilles::PSBuilder(nlep, nhad).Beam(beam, masses, 1)
                                                   .Hadron(model -> PhaseSpace(), masses)
                                                   .FinalState(T::Name(), masses, gauge_boson_mass).build();
    achilles::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
    return channel;
}

#ifdef ENABLE_BSM
template<typename T>
achilles::Channel<achilles::FourVector> BuildChannelSherpa(achilles::NuclearModel *model, size_t nlep, size_t nhad,
                                                       std::shared_ptr<achilles::Beam> beam,
                                                       const std::vector<double> &masses) {
    achilles::Channel<achilles::FourVector> channel;
    auto massesGeV = masses;
    for(auto &mass : massesGeV) mass /= (1000*1000);
    channel.mapping = achilles::PSBuilder(nlep, nhad).Beam(beam, masses, 1)
                                                   .Hadron(model -> PhaseSpace(), masses)
                                                   .SherpaFinalState(T::Name(), massesGeV).build();
    achilles::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
    return channel;
}

achilles::Channel<achilles::FourVector> BuildGenChannel(achilles::NuclearModel *model, size_t nlep, size_t nhad,
                                                    std::shared_ptr<achilles::Beam> beam,
                                                    std::unique_ptr<PHASIC::Channels> final_state,
                                                    const std::vector<double> &masses) {
    achilles::Channel<achilles::FourVector> channel;
    channel.mapping = achilles::PSBuilder(nlep, nhad).Beam(beam, masses, 1)
                                                   .Hadron(model -> PhaseSpace(), masses)
                                                   .GenFinalState(std::move(final_state)).build();
    achilles::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
    return channel;
}
#endif

achilles::EventGen::EventGen(const std::string &configFile,
                             std::vector<std::string> shargs) {
    config = YAML::LoadFile(configFile);

    // Setup random number generator
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["Seed"])
        if(config["Initialize"]["Seed"].as<int>() > 0)
            seed = config["Initialize"]["Seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Setup unweighter
    unweighter = UnweighterFactory::Initialize(config["Unweighting"]["Name"].as<std::string>(),
                                               config["Unweighting"]);

    // Load initial state, massess
    spdlog::trace("Initializing the beams");
    beam = std::make_shared<Beam>(config["Beams"].as<Beam>());
    nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());

    // Set potential for the nucleus
    auto potential_name = config["Nucleus"]["Potential"]["Name"].as<std::string>();
    auto potential = achilles::PotentialFactory::Initialize(potential_name,
                                                            config["Nucleus"]["Potential"]);
    nucleus -> SetPotential(std::move(potential));
    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Run"].as<bool>());
    if(config["Cascade"]["Run"].as<bool>()) {
        cascade = std::make_unique<Cascade>(config["Cascade"].as<Cascade>());
    } else {
        cascade = nullptr;
    }

    // Initialize the lepton final states
    spdlog::debug("Initializing the leptonic final states");
    auto leptonicProcess = config["Process"].as<achilles::Process_Info>();
    // TODO: Handle the beam initial state better
    if(beam -> BeamIDs().size() > 1)
        throw std::runtime_error("Multiple processes are not implemented yet. Please use only one beam.");
    leptonicProcess.ids.insert(leptonicProcess.ids.begin(),
                               beam -> BeamIDs().begin(), beam -> BeamIDs().end());

#ifdef ENABLE_BSM
    // Initialize sherpa processes
    achilles::SherpaMEs *sherpa = new achilles::SherpaMEs();
    std::string model = config["Process"]["Model"].as<std::string>();
    std::string param_card = config["Process"]["ParamCard"].as<std::string>();
    //if(model == "SM") model = "SM_Nuc";
    shargs.push_back("MODEL=" + model);
    shargs.push_back("UFO_PARAM_CARD=" + param_card);
    sherpa -> Initialize(shargs);
#endif

    // Initialize the nuclear model
    spdlog::debug("Initializing nuclear model");
    const auto model_name = config["NuclearModel"]["Model"].as<std::string>();
    auto nuclear_model = NuclearModelFactory::Initialize(model_name, config, nucleus);
    auto proc_group = nuclear_model -> AllowedStates(leptonicProcess);
    spdlog::debug("Process Group: {}", proc_group);

#ifdef ENABLE_BSM
    spdlog::debug("Initializing leptonic currents");
    for(size_t i = 0; i < proc_group.Processes().size(); ++i) {
        auto &process = proc_group.Process(i);
        if(!sherpa -> InitializeProcess(process)) {

            spdlog::error("Cannot initialize hard process");
            exit(1);
        }
        process.m_mom_map = sherpa -> MomentumMap(process.Ids());
    }
#else
    // Dummy call to remove unused error
    if(shargs.size() != 0)
        spdlog::warn("Sherpa arguments are being ignored since Achilles is not compiled with Sherpa support");
    shargs.clear();
    for(size_t i = 0; i < proc_group.Processes().size(); ++i) {
        auto &process = proc_group.Process(i);
        for(size_t j = 0; j < process.Ids().size(); ++j) {
            process.m_mom_map[j] = process.Ids()[j];
            spdlog::debug("PIDs are {}", process.m_mom_map[j]);
        }
    }
#endif

    // Initialize hard cross-sections
    spdlog::debug("Initializing hard interaction");
    scattering = std::make_shared<HardScattering>();
    scattering -> SetProcessGroup(proc_group);
#ifdef ENABLE_BSM
    scattering -> SetSherpa(sherpa);
#endif
    scattering -> SetNuclear(std::move(nuclear_model));

    // Setup channels
    spdlog::debug("Initializing phase space");
    // TODO: Take into account all processes in a process group
    std::vector<double> masses = scattering -> ProcessGroup().Processes()[0].Masses();
    spdlog::debug("Masses = [{}]", fmt::join(masses.begin(), masses.end(), ", "));

    if(config["TestingPS"]) {
        Channel<FourVector> channel = BuildChannelTest(config["TestingPS"], beam);
        integrand.AddChannel(std::move(channel));
    } else {
#ifndef ENABLE_BSM

        if(scattering -> ProcessGroup().Multiplicity() == 4) {
            Channel<FourVector> channel0 = BuildChannel<TwoBodyMapper>(scattering -> Nuclear(), 2, 2,
                                                                       beam, masses);
            integrand.AddChannel(std::move(channel0));
        } 

        else if(scattering -> ProcessGroup().Multiplicity() == 5) {
            double gauge_boson_mass = 80.354;
            Channel<FourVector> channel1 = BuildChannel<ThreeBodyMapper>(scattering -> Nuclear(), 2, 3,
                                                                         beam, masses, gauge_boson_mass);
            integrand.AddChannel(std::move(channel1));
        } 

        else {
            const std::string error = fmt::format("Can only handle 2->2 and 2->3 processes without "
                                                  "BSM being enabled. "
                                                  "Got a 2->{} process",
                                                  scattering -> ProcessGroup().Multiplicity()-2);
            throw std::runtime_error(error);
        }
#else

        // TODO: Take into account all processes in a process group
    auto channels = sherpa -> GenerateChannels(scattering -> ProcessGroup().Processes()[0].Ids());

    size_t count = 0;
 

        for(auto & chan : channels) {
            Channel<FourVector> channel = BuildGenChannel(scattering -> Nuclear(), 
                                                          scattering -> ProcessGroup().Processes()[0].ids.size(), 2,
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
    if(config["Main"]["DoRotate"])
        doRotate = config["Main"]["DoRotate"].as<bool>();

    // Setup Cuts
    if(config["Main"]["HardCuts"])
        doHardCuts = config["Main"]["HardCuts"].as<bool>();
    spdlog::info("Apply hard cuts? {}", doHardCuts);
    hard_cuts = config["HardCuts"].as<achilles::CutCollection>();

    // Setup outputs
    auto output = config["Main"]["Output"];
    bool zipped = true;
    if(output["Zipped"])
        zipped = output["Zipped"].as<bool>();
    spdlog::trace("Outputing as {} format", output["Format"].as<std::string>());
    if(output["Format"].as<std::string>() == "Achilles") {
        writer = std::make_unique<AchillesWriter>(output["Name"].as<std::string>(), zipped);
#ifdef ENABLE_HEPMC3
    } else if(output["Format"].as<std::string>() == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(output["Name"].as<std::string>(), zipped);
#endif
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}",
                                      output["Format"].as<std::string>());
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
        if(config["Initialize"]["Accuracy"])
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
    runCascade = config["Cascade"]["Run"].as<bool>();
    integrator.Parameters().ncalls = config["Main"]["NEvents"].as<size_t>();
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
                       config["Main"]["NEvents"].as<size_t>());
        }
    }
    // Initialize the event, which generates the nuclear configuration
    // and initializes the beam particle for the event
    Event event(nucleus, mom, wgt);

    // Initialize the particle ids for the processes
    // TODO: Take into account all processes within the group
    const auto pids = scattering -> ProcessGroup().Processes()[0].ids;

    // Setup flux value
    event.Flux() = beam -> EvaluateFlux(pids[0], mom[1]);

    spdlog::debug("Event Phase Space:");
    size_t idx = 0;
    for(const auto &momentum : event.Momentum()) {
        spdlog::debug("\t{}: {} (M2 = {})", ++idx, momentum, momentum.M2());
    }

    // Calculate the hard cross sections and select one for initial state
    spdlog::debug("Calculating cross section");

    // Obtain the cross section for the event 
    auto xsecs = scattering -> CrossSection(event);

    // Initialize the event
    spdlog::debug("Filling the event");
    if(!scattering -> FillEvent(event, xsecs)) {
        if(outputEvents) {
            event.CalcTotalCrossSection({0});
            event.CalcWeight();
            spdlog::trace("Outputting the event");
            writer -> Write(event);
            // Update number of calls needed to ensure the number of generated events
            // is the same as that requested by the user
            integrator.Parameters().ncalls++;
        }
        return 0;
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

    event.CalcWeight();
    spdlog::trace("Weight: {}", event.Weight());

    // if((event.Momentum()[3]+event.Momentum()[4]).M() < 400) {
    //     spdlog::info("Mass issue");
    //     spdlog::drop("achilles");
    //     CreateLogger(0, 5);
    // }

    // Perform hard cuts
    if(doHardCuts) {
        spdlog::debug("Making hard cuts");
        if(!MakeCuts(event)) {
            // Short-circuit the evaluation
            // We want Vegas to adapt to avoid these points, i.e.,
            // the integrand should be interpreted as zero in this region
            if(outputEvents) {
                event.CalcTotalCrossSection({0});
                event.CalcWeight();
                writer -> Write(event);
                // Update number of calls needed to ensure the number of generated events
                // is the same as that requested by the user
                integrator.Parameters().ncalls++;
            }
            return 0;
        }
    }

    // Run the cascade if needed
    if(runCascade) {
        spdlog::trace("Hadrons:");
        idx = 0;
        for(const auto &particle : event.Hadrons()) {
            if(particle.Status() == ParticleStatus::initial_state
                && particle.ID() == PID::proton())
                spdlog::trace("\t{}: {}", idx, particle);
            ++idx;
        }
        spdlog::debug("Runnning cascade");
        cascade -> Evolve(&event);

        spdlog::trace("Hadrons (Post Cascade):");
        idx = 0;
        for(const auto &particle : event.Hadrons()) {
            spdlog::trace("\t{}: {}", ++idx, particle);
        }
    } else {
        for(auto & nucleon : event.CurrentNucleus()->Nucleons()) {
            if(nucleon.Status() == ParticleStatus::propagating) {
                nucleon.Status() = ParticleStatus::final_state;
            }
        }
    }

    // Write out events
    if(outputEvents) {
        // Rotate cuts into plane of outgoing electron before writing
        if (doRotate)
            Rotate(event);
        // Perform event-level final cuts before writing
        bool outputCurrentEvent = true;

        if(outputCurrentEvent) {
            // Keep a running total of the number of surviving events
            event.Finalize();
            if(!unweighter->AcceptEvent(event.Weight())) {
                // Update number of calls needed to ensure the number of generated events
                // is the same as that requested by the user
                integrator.Parameters().ncalls++;
            }
            writer -> Write(event);
        }
    } else {
        unweighter->AddEvent(event.Weight());
    }

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
