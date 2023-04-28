#include "Achilles/EventGen.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Channels.hh"
#include "Achilles/ComplexFmt.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/HardScattering.hh"
#include "Achilles/HardScatteringFactory.hh"
#include "Achilles/Logging.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Units.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/SherpaInterface.hh"
#else
namespace achilles {
class SherpaInterface {};
} // namespace achilles
#endif

#ifdef ACHILLES_ENABLE_HEPMC3
#include "plugins/HepMC3/HepMC3EventWriter.hh"
#endif

#include "yaml-cpp/yaml.h"

achilles::EventGen::EventGen(const std::string &configFile, std::vector<std::string> shargs) {
    config = YAML::LoadFile(configFile);

    // Setup random number generator
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["Seed"])
        if(config["Initialize"]["Seed"].as<int>() > 0)
            seed = config["Initialize"]["Seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load initial beam and nucleus
    // TODO: Allow for multiple nuclei
    spdlog::trace("Initializing the beams");
    beam = std::make_shared<Beam>(config["Beams"].as<Beam>());
    nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());

    // Set potential for the nucleus
    auto potential_name = config["Nucleus"]["Potential"]["Name"].as<std::string>();
    auto potential = achilles::PotentialFactory::Initialize(potential_name, nucleus,
                                                            config["Nucleus"]["Potential"]);
    nucleus->SetPotential(std::move(potential));

    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Run"].as<bool>());
    if(config["Cascade"]["Run"].as<bool>()) {
        cascade = std::make_unique<Cascade>(config["Cascade"].as<Cascade>());
    } else {
        cascade = nullptr;
    }

    // Initialize the Nuclear models
    spdlog::debug("Initializing nuclear models");
    auto models = LoadModels(config);

    // Initialize Sherpa
    p_sherpa = new achilles::SherpaInterface();
#ifdef ACHILLES_SHERPA_INTERFACE
    auto model_name = config["SherpaOptions"]["Model"].as<std::string>();
    auto param_card = config["SherpaOptions"]["ParamCard"].as<std::string>();
    int qed = 0;
    if(config["SherpaOptions"]["QEDShower"])
        if(config["SherpaOptions"]["QEDShower"].as<bool>()) qed = 1;
    shargs.push_back(fmt::format("ME_QED={}", qed));
    if(model_name == "SM") model_name = "SM_Nuc";
    shargs.push_back("MODEL=" + model_name);
    shargs.push_back("UFO_PARAM_CARD=" + param_card);
    p_sherpa->Initialize(shargs);
#endif

    // Initialize the Processes
    spdlog::debug("Initializing processes");
    if(beam->BeamIDs().size() > 1)
        throw std::runtime_error("Multiple processes are not implemented yet. "
                                 "Please use only one beam.");
    for(auto &model : models) {
        auto groups =
            ProcessGroup::ConstructProcessGroups(config, model.second.get(), beam, nucleus);
        std::cout << groups.size() << std::endl;
        for(auto &group : groups) {
            for(const auto &process : group.second.Processes())
                spdlog::info("Found Process: {}", process.Info());
            group.second.SetupBackend(config, std::move(model.second), p_sherpa);
            process_groups.push_back(std::move(group.second));
        }
    }

    // Setup Multichannel integrators
    for(auto &group : process_groups) { group.SetupIntegration(config); }

    // Decide whether to rotate events to be measured w.r.t. the lepton plane
    if(config["Main"]["DoRotate"]) doRotate = config["Main"]["DoRotate"].as<bool>();

    // Setup Cuts
    if(config["Main"]["HardCuts"]) doHardCuts = config["Main"]["HardCuts"].as<bool>();
    spdlog::info("Apply hard cuts? {}", doHardCuts);
    hard_cuts = config["HardCuts"].as<achilles::CutCollection>();

    // Setup outputs
    auto output = config["Main"]["Output"];
    bool zipped = true;
    if(output["Zipped"]) zipped = output["Zipped"].as<bool>();
    spdlog::trace("Outputing as {} format", output["Format"].as<std::string>());
    if(output["Format"].as<std::string>() == "Achilles") {
        writer = std::make_unique<AchillesWriter>(output["Name"].as<std::string>(), zipped);
#ifdef ACHILLES_ENABLE_HEPMC3
    } else if(output["Format"].as<std::string>() == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(output["Name"].as<std::string>(), zipped);
#endif
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}",
                                      output["Format"].as<std::string>());
        throw std::runtime_error(msg);
    }
    writer->WriteHeader(configFile);
}

void achilles::EventGen::Initialize() {
    // TODO: Save trained integrators
    spdlog::info("Starting optimization runs");
    for(auto &process : process_groups) { process.Optimize(); }
}

void achilles::EventGen::GenerateEvents() {
    outputEvents = true;
    runCascade = config["Cascade"]["Run"].as<bool>();
    integrator.Parameters().ncalls = config["Main"]["NEvents"].as<size_t>();
    std::string filename = fmt::format(
        "events_{}.txt", config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>());
    outputfile.open(filename);
    integrator(integrand);
    fmt::print("\n");
    auto result = integrator.Summary();
    fmt::print("Integral = {:^8.5e} +/- {:^8.5e} ({:^8.5e} %)\n", result.results.back().Mean(),
               result.results.back().Error(),
               result.results.back().Error() / result.results.back().Mean() * 100);
    fmt::print("Unweighting efficiency: {:^8.5e} %\n", unweighter->Efficiency() * 100);
}

double achilles::EventGen::GenerateEvent(const std::vector<FourVector> &mom, const double &wgt) {
    if(outputEvents) {
        static constexpr size_t statusUpdate = 1000;
        if(unweighter->Accepted() % statusUpdate == 0) {
            fmt::print("Generated {} / {} events\r", unweighter->Accepted(),
                       config["Main"]["NEvents"].as<size_t>());
        }
    }
    // Initialize the event, which generates the nuclear configuration
    // and initializes the beam particle for the event
    Event event(nucleus, mom, wgt);

    // Initialize the particle ids for the processes
    auto pids = scattering->Process().m_leptonic.second;
    pids.insert(pids.begin(), scattering->Process().m_leptonic.first);

    // Setup flux value
    event.Flux() = beam->EvaluateFlux(pids[0], mom[1]);

    spdlog::debug("Event Phase Space:");
    size_t idx = 0;
    for(const auto &momentum : event.Momentum()) {
        spdlog::debug("\t{}: {} (M2 = {})", ++idx, momentum, momentum.M2());
    }

    // Calculate the hard cross sections and select one for initial state
    spdlog::debug("Calculating cross section");

    // Obtain the cross section for the event
    auto xsecs = scattering->CrossSection(event);

    // Initialize the event
    spdlog::debug("Filling the event");
    if(!scattering->FillEvent(event, xsecs)) {
        if(outputEvents) {
            event.SetMEWeight(0);
            event.CalcWeight();
            spdlog::trace("Outputting the event");
            writer->Write(event);
            // Update number of calls needed to ensure the number of generated
            // events is the same as that requested by the user
            integrator.Parameters().ncalls++;
        }
        return 0;
    }

    spdlog::trace("Leptons:");
    idx = 0;
    for(const auto &particle : event.Leptons()) { spdlog::trace("\t{}: {}", ++idx, particle); }

    spdlog::trace("Hadrons:");
    idx = 0;
    for(const auto &particle : event.Hadrons()) { spdlog::trace("\t{}: {}", ++idx, particle); }

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
                event.SetMEWeight(0);
                event.CalcWeight();
                writer->Write(event);
                // Update number of calls needed to ensure the number of
                // generated events is the same as that requested by the user
                integrator.Parameters().ncalls++;
            }
            return 0;
        }
    }

    // TODO: Move to after unweighting?
    // Run the cascade if needed
    if(runCascade) {
        spdlog::trace("Hadrons:");
        idx = 0;
        for(const auto &particle : event.Hadrons()) {
            if(particle.Status() == ParticleStatus::initial_state && particle.ID() == PID::proton())
                spdlog::trace("\t{}: {}", idx, particle);
            ++idx;
        }
        spdlog::debug("Runnning cascade");
        cascade->Evolve(&event);

        spdlog::trace("Hadrons (Post Cascade):");
        idx = 0;
        for(const auto &particle : event.Hadrons()) { spdlog::trace("\t{}: {}", ++idx, particle); }
    } else {
        for(auto &nucleon : event.CurrentNucleus()->Nucleons()) {
            if(nucleon.Status() == ParticleStatus::propagating) {
                nucleon.Status() = ParticleStatus::final_state;
            }
        }
    }

    // Write out events
    if(outputEvents) {
#ifdef ACHILLES_SHERPA_INTERFACE
        // Running Sherpa interface if requested
        // Only needed when generating events and not optimizing the
        // multichannel
        // TODO: Move to after unweighting?
        if(runDecays && event.Weight() > 0) { p_sherpa->GenerateEvent(event); }
#endif

        // Rotate cuts into plane of outgoing electron before writing
        if(doRotate) Rotate(event);
        // Perform event-level final cuts before writing
        bool outputCurrentEvent = true;
        // if(doEventCuts){
        //     spdlog::debug("Making event cuts");
        //     outputCurrentEvent = MakeEventCuts(event);
        // }

        if(outputCurrentEvent) {
            // Keep a running total of the number of surviving events
            event.Finalize();
            if(!unweighter->AcceptEvent(event.Weight())) {
                // Update number of calls needed to ensure the number of
                // generated events is the same as that requested by the user
                integrator.Parameters().ncalls++;
            }
            writer->Write(event);
        } else {
            unweighter->AddEvent(event.Weight());
        }
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
    for(const auto &particle : event.Particles()) {
        if(particle.ID() == PID::electron() && particle.IsFinal()) {
            phi = particle.Momentum().Phi();
        }
    }
    // Rotate the coordiantes of particles so that all azimuthal angles phi are
    // measured with respect to the leptonic plane
    std::array<double, 9> rotation = {cos(phi), sin(phi), 0, -sin(phi), cos(phi), 0, 0, 0, 1};
    event.Rotate(rotation);
}
