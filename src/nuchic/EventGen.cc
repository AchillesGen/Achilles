#include "nuchic/EventGen.hh"
#include "nuchic/Event.hh"
#include "nuchic/EventWriter.hh"
#include "nuchic/HardScatteringFactory.hh"
#include "nuchic/HardScattering.hh"
#include "nuchic/Logging.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Beams.hh"
#include "nuchic/Cascade.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Units.hh"
#include "nuchic/ProcessInfo.hh"
#include "nuchic/NuclearModel.hh"
#include "nuchic/ComplexFmt.hh"
#include "nuchic/Units.hh"

// TODO: Turn this into a factory to reduce the number of includes
#include "nuchic/PhaseSpaceBuilder.hh"
#include "nuchic/BeamMapper.hh"
#include "nuchic/HadronicMapper.hh"
#include "nuchic/FinalStateMapper.hh"
#include "nuchic/PhaseSpaceMapper.hh"
#include "nuchic/QuasielasticTestMapper.hh"

#ifdef ENABLE_BSM
#include "plugins/Sherpa/Channels1.hh"
#include "plugins/Sherpa/Channels3.hh"
#include "plugins/Sherpa/SherpaMEs.hh"
#endif

#ifdef ENABLE_HEPMC3
#include "plugins/HepMC3/HepMC3EventWriter.hh"
#endif

#include "yaml-cpp/yaml.h"

nuchic::Channel<nuchic::FourVector> BuildChannelTest(const YAML::Node &node, std::shared_ptr<nuchic::Beam> beam) {
    nuchic::Channel<nuchic::FourVector> channel;
    channel.mapping = std::make_unique<nuchic::QuasielasticTestMapper>(node, beam);
    nuchic::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = nuchic::Vegas(map, {});
    return channel;
}

template<typename T>
nuchic::Channel<nuchic::FourVector> BuildChannel(nuchic::NuclearModel *model, size_t nlep, size_t nhad,
                                                 std::shared_ptr<nuchic::Beam> beam,
                                                 const std::vector<double> &masses) {
    nuchic::Channel<nuchic::FourVector> channel;
    channel.mapping = nuchic::PSBuilder(nlep, nhad).Beam(beam, 1)
                                                   .Hadron(model -> PhaseSpace(), masses)
                                                   .FinalState(T::Name(), masses).build();
    nuchic::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = nuchic::Vegas(map, nuchic::VegasParams{});
    return channel;
}

#ifdef ENABLE_BSM
template<typename T>
nuchic::Channel<nuchic::FourVector> BuildChannelSherpa(nuchic::NuclearModel *model, size_t nlep, size_t nhad,
                                                       std::shared_ptr<nuchic::Beam> beam,
                                                       const std::vector<double> &masses) {
    nuchic::Channel<nuchic::FourVector> channel;
    auto massesGeV = masses;
    for(auto &mass : massesGeV) mass /= (1000*1000);
    channel.mapping = nuchic::PSBuilder(nlep, nhad).Beam(beam, 1)
                                                   .Hadron(model -> PhaseSpace(), masses)
                                                   .SherpaFinalState(T::Name(), massesGeV).build();
    nuchic::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = nuchic::Vegas(map, nuchic::VegasParams{});
    return channel;
}

nuchic::Channel<nuchic::FourVector> BuildGenChannel(nuchic::NuclearModel *model, size_t nlep, size_t nhad,
                                                    std::shared_ptr<nuchic::Beam> beam,
                                                    std::unique_ptr<PHASIC::Channels> final_state,
                                                    const std::vector<double> &masses) {
    nuchic::Channel<nuchic::FourVector> channel;
    channel.mapping = nuchic::PSBuilder(nlep, nhad).Beam(beam, 1)
                                                   .Hadron(model -> PhaseSpace(), masses)
                                                   .GenFinalState(std::move(final_state)).build();
    nuchic::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = nuchic::Vegas(map, nuchic::VegasParams{});
    return channel;
}
#endif

nuchic::EventGen::EventGen(const std::string &configFile, SherpaMEs *const sherpa) :
  runCascade{false}, outputEvents{false} {
    config = YAML::LoadFile(configFile);

    // Setup random number generator
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["Seed"])
        if(config["Initialize"]["Seed"].as<int>() > 0)
            seed = config["Initialize"]["Seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load initial state, massess
    beam = std::make_shared<Beam>(config["Beams"].as<Beam>());
    nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());

    // Set potential for the nucleus
    auto potential_name = config["Nucleus"]["Potential"]["Name"].as<std::string>();
    auto potential = nuchic::PotentialFactory::Initialize(potential_name,
                                                          nucleus,
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
    auto leptonicProcess = config["Process"].as<nuchic::Process_Info>();
    // TODO: Handle the beam initial state better
    if(beam -> BeamIDs().size() > 1)
        throw std::runtime_error("Multiple processes are not implemented yet. Please use only one beam.");
    leptonicProcess.m_ids.insert(leptonicProcess.m_ids.begin(),
                                 beam -> BeamIDs().begin(), beam -> BeamIDs().end());

    // Initialize the nuclear model
    spdlog::debug("Initializing nuclear model");
    const auto model_name = config["NuclearModel"]["Model"].as<std::string>();
    auto nuclear_model = NuclearModelFactory::Initialize(model_name, config);
    nuclear_model -> AllowedStates(leptonicProcess);
    spdlog::debug("Process: {}", leptonicProcess);

    // Initialize sherpa processes
    spdlog::debug("Initializing leptonic currents");
    if(!sherpa->InitializeProcess(leptonicProcess)) {
        spdlog::error("Cannot initialize hard process");
        exit(1);
    }
    leptonicProcess.m_mom_map = sherpa -> MomentumMap(leptonicProcess.Ids());

    // Initialize hard cross-sections
    spdlog::debug("Initializing hard interaction");
    scattering = std::make_shared<HardScattering>();
    scattering -> SetProcess(leptonicProcess);
    scattering -> SetSherpa(sherpa);
    scattering -> SetNuclear(std::move(nuclear_model));

    // Setup channels
    spdlog::debug("Initializing phase space");
    std::vector<double> masses = scattering -> Process().Masses();
    spdlog::trace("Masses = [{}]", fmt::join(masses.begin(), masses.end(), ", "));
    if(config["TestingPS"]) {
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
        auto channels = sherpa -> GenerateChannels(scattering -> Process().Ids());
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
    if(config["Main"]["DoRotate"])
        doRotate = config["Main"]["DoRotate"].as<bool>();

    // Setup Cuts
    if(config["Main"]["HardCuts"])
        doHardCuts = config["Main"]["HardCuts"].as<bool>();
    spdlog::info("Apply hard cuts? {}", doHardCuts);
    hard_cuts = config["HardCuts"].as<nuchic::CutCollection>();

    // if(config["Main"]["EventCuts"])
    //     doEventCuts = config["Main"]["EventCuts"].as<bool>();
    // spdlog::info("Apply event cuts? {}", doEventCuts);
    // event_cuts = config["EventCuts"].as<nuchic::CutCollection>();

    // Setup outputs
    auto output = config["Main"]["Output"];
    bool zipped = true;
    if(output["Zipped"])
        zipped = output["Zipped"].as<bool>();
    if(output["Format"].as<std::string>() == "Nuchic") {
        writer = std::make_unique<NuchicWriter>(output["Name"].as<std::string>(), zipped);
#ifdef ENABLE_HEPMC3
    } else if(output["Format"].as<std::string>() == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(output["Name"].as<std::string>(), zipped);
#endif
    }
    writer -> WriteHeader(configFile);

    hist = Histogram(200, 0.0, 400.0, "energy");
    hist2 = Histogram(100, 0.0, 800.0, "momentum");
    hist3 = Histogram(50, -1.0, 1.0, "angle");
    hist4 = Histogram(200, 320.0, 520.0, "invariant_mass");
    hist5 = Histogram(300, 0.0, 300.0, "omega");
    hist6 = Histogram(200, 4.0, 8.0, "wgt");
}

void nuchic::EventGen::Initialize() {
    // TODO: Clean up loading of previous results
    auto func = [&](const std::vector<FourVector> &mom, const double &wgt) {
        return GenerateEvent(mom, wgt);
    };
    try {
        YAML::Node old_results = YAML::LoadFile("results.yml");
        integrator = old_results["Multichannel"].as<MultiChannel>();
        integrand = old_results["Channels"].as<Integrand<FourVector>>();
        YAML::Node results;
        results["Multichannel"] = integrator;
        results["Channels"] = integrand;
        integrand.Function() = func;
    } catch(const YAML::BadFile &e) {
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
    }
}

void nuchic::EventGen::GenerateEvents() {
    outputEvents = true;
    runCascade = config["Cascade"]["Run"].as<bool>();
    integrator.Parameters().ncalls = config["Main"]["NEvents"].as<size_t>();
    integrator(integrand);
    fmt::print("\n");
    auto result = integrator.Summary();
    fmt::print("Integral = {:^8.5e} +/- {:^8.5e} ({:^8.5e} %)\n",
               result.results.back().Mean(), result.results.back().Error(),
               result.results.back().Error() / result.results.back().Mean()*100);

    hist.Save(config["HistTest1"].as<std::string>());
    hist2.Save(config["HistTest2"].as<std::string>());
    hist3.Save(config["HistTest3"].as<std::string>());
    hist4.Save(config["HistTest4"].as<std::string>());
    hist5.Save(config["HistTest5"].as<std::string>());
    hist6.Save(config["HistTest6"].as<std::string>());
}

double nuchic::EventGen::GenerateEvent(const std::vector<FourVector> &mom, const double &wgt) {
    if(outputEvents) {
        static size_t ievent = 0;
        constexpr size_t statusUpdate = 10000;
        if(++ievent % statusUpdate == 0) {
            fmt::print("Generated {} / {} events\r", ievent, config["Main"]["NEvents"].as<size_t>());
        }
    }
    // Initialize the event, which generates the nuclear configuration
    // and initializes the beam particle for the event
    Event event(nucleus, mom, wgt);

    // Initialize the particle ids for the processes
    const auto pids = scattering -> Process().m_ids;

    // Calculate the hard cross sections and select one for initial state
    spdlog::debug("Calculating cross section");

    // Obtain the cross section for the event 
    auto xsecs = scattering -> CrossSection(event);

    // Initialize the event
    if(!scattering -> FillEvent(event, xsecs)) {
        if(outputEvents) {
            event.SetMEWeight(0);
            spdlog::trace("Outputting the event");
            writer -> Write(event);
        }
        return 0;
    }

    spdlog::trace("Event Phase Space:");
    size_t idx = 0;
    for(const auto &momentum : event.Momentum()) {
        spdlog::trace("\t{}: {}", ++idx, momentum);
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

    spdlog::trace("Weight: {}", event.Weight());

    // if((event.Momentum()[3]+event.Momentum()[4]).M() < 400) {
    //     spdlog::info("Mass issue");
    //     spdlog::drop("nuchic");
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
                event.SetMEWeight(0);
                writer -> Write(event);
            }
            return 0;
        }
    }

    // Run the cascade if needed
    if(runCascade) {
        spdlog::debug("Runnning cascade");
        cascade -> Evolve(&event);
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
        // if(doEventCuts){
        //     spdlog::debug("Making event cuts");
        //     outputCurrentEvent = MakeEventCuts(event);
        // }

        if(outputCurrentEvent) {
            // Keep a running total of the number of surviving events
            event.Finalize();
            writer -> Write(event);
            const auto energy = (Constant::mN - event.Momentum()[0].E());
            const auto calls = static_cast<double>(integrator.Parameters().ncalls);
            hist.Fill(energy, event.Weight()/calls);
            const auto momentum = event.Momentum()[0].P();
            hist2.Fill(momentum, event.Weight()/calls);
            const auto cosTheta = event.Momentum()[2].CosTheta();
            hist3.Fill(cosTheta, event.Weight()/calls);
            const auto energy_lepton = (event.Momentum()[3]+event.Momentum()[4]).M();
            hist4.Fill(energy_lepton, event.Weight()/calls);
            const auto omega = event.Momentum()[1].E() - (event.Momentum()[3]+event.Momentum()[4]).E();
            hist5.Fill(omega, event.Weight()/calls/(2*M_PI));
            hist6.Fill(log10(event.Weight()));
        }
    }

    // Always return the weight when the event passes the initial hard cut.
    // Even if events do not survive the final event-level cuts, Vegas should
    // still interpret the integrand as nonzero in this region.
    return event.Weight();
}

bool nuchic::EventGen::MakeCuts(Event &event) {
    return hard_cuts.EvaluateCuts(event.Particles());
}

// TODO: Create Analysis level cuts
/*
bool nuchic::EventGen::MakeEventCuts(Event &event) {
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

void nuchic::EventGen::Rotate(Event &event) {
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
