#include "Achilles/EventGen.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Channels.hh"
#include "Achilles/Debug.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventWriter.hh"
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
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdouble-promotion"
#endif
#include "plugins/HepMC3/HepMC3EventWriter.hh"
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#include "plugins/NuHepMC/NuHepMCWriter.hh"
#endif

#include "fmt/std.h"
#include "yaml-cpp/yaml.h"

achilles::EventGen::EventGen(const std::string &configFile, std::vector<std::string> shargs)
    : config{configFile} {
    // Turning off decays in Sherpa. This is a temporary fix until we can get ISR and FSR properly
    // working in SHERPA.
    runDecays = false;

    // Setup random number generator
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config.Exists("Options/Initialize/Seed"))
        if(config.GetAs<int>("Options/Initialize/Seed") > 0)
            seed = config.GetAs<unsigned int>("Options/Initialize/Seed");
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load initial beam and nuclei
    // TODO: Allow for multiple nuclei
    spdlog::trace("Initializing the beams");
    beam = std::make_shared<Beam>(config.GetAs<Beam>("Beams"));
    nuclei = config.GetAs<std::vector<std::shared_ptr<Nucleus>>>("Nuclei");

    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config.GetAs<bool>("Cascade/Run"));
    if(config.GetAs<bool>("Cascade/Run")) {
        cascade = std::make_unique<Cascade>(config.GetAs<Cascade>("Cascade"));
    } else {
        cascade = nullptr;
    }

    // Initialize decays
    runDecays = config.GetAs<bool>("Main/RunDecays");

    // Initialize Sherpa
    auto backend = config.GetAs<std::string>("Backend/Name");
    if(backend.find("Sherpa") != std::string::npos || backend.find("BSM") != std::string::npos) {
#ifdef ACHILLES_SHERPA_INTERFACE
        p_sherpa = new achilles::SherpaInterface();
        std::string model_name = "SMnu";
        if(config.Exists("SherpaOptions/Model"))
            model_name = config.GetAs<std::string>("SherpaOptions/Model");
        auto param_card = config.GetAs<std::string>("SherpaOptions/ParamCard");
        int qed = 0;
        if(config.Exists("SherpaOptions/QEDShower"))
            if(config.GetAs<bool>("SherpaOptions/QEDShower")) qed = 1;
        shargs.push_back(fmt::format("ME_QED={}", qed));
        if(model_name == "SM") model_name = "SMnu";
        shargs.push_back("MODEL=" + model_name);
        shargs.push_back("UFO_PARAM_CARD=" + param_card);
        p_sherpa->Initialize(shargs);
#else
        shargs.clear();
        throw std::runtime_error("Achilles has not been compiled with Sherpa support!");
#endif
    } else
        p_sherpa = nullptr;

    // Initialize the Processes
    spdlog::debug("Initializing processes");
    // if(beam->BeamIDs().size() > 1)
    //     throw std::runtime_error("Multiple processes are not implemented yet. "
    //                              "Please use only one beam.");

    for(const auto &nucleus : nuclei) {
        // Initialize the Nuclear models for each nuclei
        spdlog::debug("Initializing nuclear models");
        // TODO: Convert to using Settings object instead of calling Root()
        auto models = LoadModels(config.Root());
        for(auto &model : models) {
            auto groups =
                ProcessGroup::ConstructGroups(config.Root(), model.second.get(), beam, nucleus);
            for(auto &group : groups) {
                for(const auto &process : group.second.Processes())
                    spdlog::info("Found Process: {}", process.Info());
                group.second.SetupBackend(config.Root(), std::move(model.second), p_sherpa);
                process_groups.push_back(std::move(group.second));
            }
        }
    }

    // Setup Multichannel integrators and remove invalid configurations
    // TODO: Convert to using Settings object instead of calling Root()
    process_groups.erase(
        std::remove_if(process_groups.begin(), process_groups.end(),
                       [&](ProcessGroup &group) { return !group.SetupIntegration(config.Root()); }),
        process_groups.end());

    for(auto &group : process_groups) {
        std::cout << std::hex << std::hash<ProcessGroup>{}(group) << std::endl;
    }

    // Decide whether to rotate events to be measured w.r.t. the lepton plane
    if(config.Exists("Main/DoRotate")) doRotate = config.GetAs<bool>("Main/DoRotate");

    // Setup Cuts
    if(config.Exists("Main/HardCuts")) {
        doHardCuts = config.GetAs<bool>("Main/HardCuts");
        spdlog::info("Apply hard cuts? {}", doHardCuts);
        hard_cuts = config.GetAs<achilles::CutCollection>("HardCuts");
        for(auto &group : process_groups) group.SetCuts(hard_cuts);
    }

    // Setup outputs
    bool zipped = true;
    if(config.Exists("Main/Output/Zipped")) zipped = config.GetAs<bool>("Main/Output/Zipped");
    auto format = config.GetAs<std::string>("Main/Output/Format");
    auto name = config.GetAs<std::string>("Main/Output/Name");
    spdlog::trace("Outputing as {} format", format);
    if(format == "Achilles") {
        writer = std::make_unique<AchillesWriter>(name, zipped);
#ifdef ACHILLES_ENABLE_HEPMC3
    } else if(format == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(name, zipped);
    } else if(format == "NuHepMC") {
        writer = std::make_unique<NuHepMCWriter>(name, zipped);
#endif
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}", format);
        throw std::runtime_error(msg);
    }
    writer->WriteHeader(configFile, process_groups);
}

void achilles::EventGen::Initialize() {
    if(config.Exists("Backend/Options/DebugEvents")) {
        auto &group = process_groups[0];
        DebugEvents events(config.GetAs<std::string>("Backend/Options/DebugEvents"),
                           group.Multiplicity());
        for(const auto &evt : events.events) { group.SingleEvent(evt, 1); }
        exit(0);
    }

    // Initialize cache
    Filesystem::Cache cache;
    if(config.Exists("Cache/Path")) cache.Path() = config.GetAs<std::string>("Cache/Path");

    if(!config.Exists("Cache/Load") || config.GetAs<bool>("Cache/Load")) {
        for(auto &group : process_groups) {
            if(cache.FindCachedState(std::hash<ProcessGroup>{}(group))) {
                // spdlog::debug("Loading cached state for {}", group);
                if(!cache.LoadState(group)) {
                    // TODO: Figure out how to stringify process group
                    spdlog::warn("Failed to load cached state for {}",
                                 std::hash<ProcessGroup>{}(group));
                }
            }
        }
    }

    spdlog::info("Starting optimization runs");
    for(auto &group : process_groups) {
        group.Optimize();
        m_group_weights.push_back(group.MaxWeight());
        m_max_weight += group.MaxWeight();
        spdlog::info("Group weights: {} / {}",
                     fmt::join(m_group_weights.begin(), m_group_weights.end(), ", "), m_max_weight);
        std::cout << "Estimated unweighting eff for this group: ";
        for(auto &process : group.Processes()) std::cout << process.UnweightEff() << " ";
        std::cout << std::endl;

        if(!config.Exists("Cache/Save") || config.GetAs<bool>("Cache/Save")) {
            if(!cache.SaveState(group)) {
                // TODO: Figure out how to stringify process group
                spdlog::warn("Failed to save cached state for {}",
                             std::hash<ProcessGroup>{}(group));
            }
        }
    }
    for(auto &wgt : m_group_weights) wgt /= m_max_weight;
    spdlog::info("Group weights: {}",
                 fmt::join(m_group_weights.begin(), m_group_weights.end(), ", "));
    spdlog::info("Finished optimization");
}

void achilles::EventGen::GenerateEvents() {
    outputEvents = true;
    runCascade = cascade != nullptr;

    const auto nevents = config["Main/NEvents"].as<size_t>();
    size_t accepted = 0;
    while(accepted < nevents) {
        static constexpr size_t statusUpdate = 1000;
        if(accepted % statusUpdate == 0) {
            fmt::print("Generated {} / {} events\r", accepted, nevents);
        }
        if(GenerateSingleEvent()) accepted++;
    }
}

bool achilles::EventGen::GenerateSingleEvent() {
    // Select the process group and generate an event
    auto &group = process_groups[Random::Instance().SelectIndex(m_group_weights)];
    auto &&event = group.GenerateEvent();
    event.Weight() *= m_max_weight;
    if(event.Weight() == 0) {
        writer->Write(event);
        return false;
    }
    if(spdlog::get_level() == spdlog::level::trace) event.Display();

    auto init_nuc = group.GetNucleus()->InitParticle();
    std::vector<Particle> init_parts;
    for(const auto &nucleon : event.Hadrons()) {
        if(nucleon.Status() == ParticleStatus::initial_state) { init_parts.push_back(nucleon); }
    }
    // TODO: Handle multiple positions from MEC
    event.History().AddVertex(init_parts[0].Position(), {init_nuc}, init_parts,
                              EventHistory::StatusCode::target);
    // Setup beam in history
    auto init_lep = event.Leptons()[0];
    auto init_beam = init_lep;
    init_beam.Status() = ParticleStatus::beam;
    const double max_energy = beam->MaxEnergy();
    init_beam.Momentum() = {max_energy, 0, 0, max_energy};
    event.History().AddVertex({}, {init_beam}, {init_lep}, EventHistory::StatusCode::beam);

    // TODO: Figure out how to best handle tracking this with the cascade and decays
    std::vector<Particle> primary_out, propagating;
    for(const auto &part : event.Particles()) {
        if(part.IsFinal()) primary_out.push_back(part);
        if(part.IsPropagating()) {
            primary_out.push_back(part);
            propagating.push_back(part);
        }
    }
    init_parts.push_back(init_lep);
    event.History().AddVertex(init_parts[0].Position(), init_parts, primary_out,
                              EventHistory::StatusCode::primary);

    // TODO: Determine if cascade or Sherpa Decays go first
    // Cascade the nucleus
    if(runCascade) {
        spdlog::debug("Runnning cascade");
        cascade->Evolve(event, group.GetNucleus());

#ifdef ACHILLES_EVENT_DETAILS
        spdlog::trace("Hadrons (Post Cascade):");
        size_t idx = 0;
        for(const auto &particle : event.Hadrons()) { spdlog::trace("\t{}: {}", ++idx, particle); }
#endif
    } else {
        for(auto &nucleon : event.Hadrons()) {
            if(nucleon.Status() == ParticleStatus::propagating) {
                nucleon.Status() = ParticleStatus::final_state;
            }
        }
    }
    // Update particle statuses in history to account for after the cascade
    event.History().UpdateStatuses(event.Hadrons());

#ifdef ACHILLES_SHERPA_INTERFACE
    // Running Sherpa interface if requested
    if(runDecays) { p_sherpa->GenerateEvent(event); }
#endif

    writer->Write(event);
    return true;
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
