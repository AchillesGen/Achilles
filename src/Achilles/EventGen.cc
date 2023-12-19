#include "Achilles/EventGen.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Channels.hh"
#include "Achilles/ComplexFmt.hh"
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
#include "plugins/HepMC3/HepMC3EventWriter.hh"
#include "plugins/NuHepMC/NuHepMCWriter.hh"
#endif

#include "yaml-cpp/yaml.h"

achilles::EventGen::EventGen(const std::string &configFile, std::vector<std::string> shargs) {
    config = YAML::LoadFile(configFile);

    // Turning off decays in Sherpa. This is a temporary fix until we can get ISR and FSR properly
    // working in SHERPA.
    runDecays = false;

    // Setup random number generator
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["Seed"])
        if(config["Initialize"]["Seed"].as<int>() > 0)
            seed = config["Initialize"]["Seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load initial beam and nuclei
    // TODO: Allow for multiple nuclei
    spdlog::trace("Initializing the beams");
    beam = std::make_shared<Beam>(config["Beams"].as<Beam>());
    nuclei = config["Nuclei"].as<std::vector<std::shared_ptr<Nucleus>>>();

    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Run"].as<bool>());
    if(config["Cascade"]["Run"].as<bool>()) {
        cascade = std::make_unique<Cascade>(config["Cascade"].as<Cascade>());
    } else {
        cascade = nullptr;
    }

    // Initialize decays
    runDecays = config["Main"]["RunDecays"].as<bool>();

    // Initialize Sherpa
    if(config["Backend"]["Name"].as<std::string>().find("Sherpa") != std::string::npos ||
       config["Backend"]["Name"].as<std::string>().find("BSM") != std::string::npos) {
#ifdef ACHILLES_SHERPA_INTERFACE
        p_sherpa = new achilles::SherpaInterface();
        std::string model_name = "SMnu";
        if(config["SherpaOptions"]["Model"])
            model_name = config["SherpaOptions"]["Model"].as<std::string>();
        auto param_card = config["SherpaOptions"]["ParamCard"].as<std::string>();
        int qed = 0;
        if(config["SherpaOptions"]["QEDShower"])
            if(config["SherpaOptions"]["QEDShower"].as<bool>()) qed = 1;
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
        auto models = LoadModels(config);
        for(auto &model : models) {
            auto groups = ProcessGroup::ConstructGroups(config, model.second.get(), beam, nucleus);
            for(auto &group : groups) {
                for(const auto &process : group.second.Processes())
                    spdlog::info("Found Process: {}", process.Info());
                group.second.SetupBackend(config, std::move(model.second), p_sherpa);
                process_groups.push_back(std::move(group.second));
            }
        }
    }

    // Setup Multichannel integrators and remove invalid configurations
    process_groups.erase(
        std::remove_if(process_groups.begin(), process_groups.end(),
                       [&](ProcessGroup &group) { return !group.SetupIntegration(config); }),
        process_groups.end());

    // Decide whether to rotate events to be measured w.r.t. the lepton plane
    if(config["Main"]["DoRotate"]) doRotate = config["Main"]["DoRotate"].as<bool>();

    // Setup Cuts
    if(config["Main"]["HardCuts"]) doHardCuts = config["Main"]["HardCuts"].as<bool>();
    spdlog::info("Apply hard cuts? {}", doHardCuts);
    hard_cuts = config["HardCuts"].as<achilles::CutCollection>();
    for(auto &group : process_groups) group.SetCuts(hard_cuts);

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
    } else if(output["Format"].as<std::string>() == "NuHepMC") {
        writer = std::make_unique<NuHepMCWriter>(output["Name"].as<std::string>(), zipped);
#endif
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}",
                                      output["Format"].as<std::string>());
        throw std::runtime_error(msg);
    }
    writer->WriteHeader(configFile, process_groups);
}

void achilles::EventGen::Initialize() {
    if(config["Backend"]["Options"]["DebugEvents"]) {
        auto &group = process_groups[0];
        DebugEvents events(config["Backend"]["Options"]["DebugEvents"].as<std::string>(),
                           group.Multiplicity());
        for(const auto &evt : events.events) { group.SingleEvent(evt, 1); }
        exit(0);
    }
    // TODO: Save trained integrators
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
    }
    for(auto &wgt : m_group_weights) wgt /= m_max_weight;
    spdlog::info("Group weights: {}",
                 fmt::join(m_group_weights.begin(), m_group_weights.end(), ", "));
    spdlog::info("Finished optimization");
}

void achilles::EventGen::GenerateEvents() {
    outputEvents = true;
    runCascade = config["Cascade"]["Run"].as<bool>();
    const auto nevents = config["Main"]["NEvents"].as<size_t>();
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
    if(event.Weight() == 0) {
        writer->Write(event);
        return false;
    }
    if(spdlog::get_level() == spdlog::level::trace) event.Display();

    auto init_nuc = event.CurrentNucleus()->InitParticle();
    std::vector<Particle> init_parts;
    for(const auto &nucleon : event.CurrentNucleus()->Nucleons()) {
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
        cascade->Evolve(&event);

        spdlog::trace("Hadrons (Post Cascade):");
        size_t idx = 0;
        for(const auto &particle : event.Hadrons()) { spdlog::trace("\t{}: {}", ++idx, particle); }
    } else {
        for(auto &nucleon : event.CurrentNucleus()->Nucleons()) {
            if(nucleon.Status() == ParticleStatus::propagating) {
                nucleon.Status() = ParticleStatus::final_state;
            }
        }
    }

#ifdef ACHILLES_SHERPA_INTERFACE
    // Running Sherpa interface if requested
    if(runDecays) { p_sherpa->GenerateEvent(event); }
#endif

    // TODO: Figure out how to best handle tracking this with the cascade and decays
    std::vector<Particle> final_part;
    for(const auto &part : event.Particles()) {
        if(part.IsFinal()) final_part.push_back(part);
    }
    init_parts.push_back(init_lep);
    event.History().AddVertex(init_parts[0].Position(), propagating, final_part,
                              EventHistory::StatusCode::cascade);

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
