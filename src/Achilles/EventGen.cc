// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/EventGen.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Channels.hh"
#include "Achilles/Debug.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/Exception.hh"
#include "Achilles/Logging.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Progress.hh"
#include "Achilles/ReferenceHandler.hh"
#include "Achilles/Units.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "Plugins/Sherpa/SherpaInterface.hh"
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
#include "Plugins/HepMC3/HepMC3EventWriter.hh"
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#include "Plugins/NuHepMC/NuHepMCWriter.hh"
#endif

#include "fmt/std.h"
#include "yaml-cpp/yaml.h"

achilles::EventGen::EventGen(const std::string &configFile, std::vector<std::string> shargs)
    : config{configFile} {
    // Turning off decays in Sherpa. This is a temporary fix until we can get ISR and FSR properly
    // working in SHERPA.
    runDecays = false;

    /** If runcard specifies Main/LogFile, output all subsequent logs to this file.
     *   Continue outputting to std::cout either way. */
    // if(config.Exists("Main/LogFile")) { (change logger) }

    // Setup random number generator
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config.Exists("Options/Initialize/Seed"))
        if(config.GetAs<int>("Options/Initialize/Seed") > 0) {
            seed = config.GetAs<unsigned int>("Options/Initialize/Seed");
            m_seed_from_card = true;
        }
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);
    m_seed = seed;

    // Load initial beam and nuclei
    spdlog::trace("Initializing the beams");
    beam = std::make_shared<Beam>(config.GetAs<Beam>("Beams"));
    nuclei = config.GetAs<std::vector<std::shared_ptr<Nucleus>>>("Nuclei");

    // Initialize Cascade parameters
    runCascade = config.GetAs<bool>("Cascade/Run");
    spdlog::debug("Cascade mode: {}", runCascade);
    cascade = runCascade ? (std::make_unique<Cascade>(config.GetAs<Cascade>("Cascade"))) : nullptr;

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
        auto models = LoadModels(config);
        for(auto &model : models) {
            auto groups = ProcessGroup::ConstructGroups(config, model.second.get(), beam, nucleus);
            for(auto &group : groups) {
                for(const auto &process : group.second.Processes())
                    spdlog::debug("Found Process: {}", process.Info());
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

    for(size_t i = 0; i < process_groups.size(); ++i)
        process_groups[i].SetId(fmt::format("G{}", i + 1));

    std::vector<ProcessMetadata> metadata;
    for(auto &group : process_groups) {
        spdlog::trace("{} has hash {:x}", group, std::hash<ProcessGroup>{}(group));
        auto group_data = group.Metadata();
        metadata.insert(metadata.end(), group_data.begin(), group_data.end());
    }

    for(const auto &data : metadata) {
        achilles::Reference ref{achilles::ReferenceType::inspire, data.inspireHEP, data.name};
        ReferenceHandler::Handle().AddReference(ref);
    }

    // Decide whether to rotate events to be measured w.r.t. the lepton plane
    if(config.Exists("Main/DoRotate")) doRotate = config.GetAs<bool>("Main/DoRotate");

    // Setup Cuts
    if(config.Exists("Main/HardCuts")) {
        doHardCuts = config.GetAs<bool>("Main/HardCuts");
        spdlog::debug("Apply hard cuts? {}", doHardCuts);
        if(doHardCuts) {
            hard_cuts = config.GetAs<achilles::CutCollection>("HardCuts");
            for(auto &group : process_groups) group.SetCuts(hard_cuts);
        }
    }

    // Setup outputs
    bool zipped =
        config.Exists("Main/Output/Zipped") ? config.GetAs<bool>("Main/Output/Zipped") : true;
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
    RunTracker::Instance().SetOutputFile(name);

    PrintConfiguration(configFile, fmt::format("{}{}", format, zipped ? ", compressed" : ""), name);
    PrintProcesses();
}

void achilles::EventGen::PrintConfiguration(const std::string &runcard, const std::string &format,
                                            const std::string &output) {
    const std::string sep = console::GetStyle().unicode ? " · " : ", ";
    std::vector<std::string> beams;
    for(const auto &pid : beam->BeamIDs()) {
        const auto flux = (*beam)[pid];
        const std::string energy =
            flux->MinEnergy() == flux->MaxEnergy()
                ? fmt::format("{:g} MeV", flux->MaxEnergy())
                : fmt::format("{:g} to {:g} MeV", flux->MinEnergy(), flux->MaxEnergy());
        beams.push_back(ParticleInfo(pid).Name() + sep + flux->Type() + sep + energy);
    }

    std::vector<std::string> targets, models;
    for(const auto &nucleus : nuclei) targets.push_back(nucleus->ToString());
    for(const auto &group : process_groups) {
        if(std::find(models.begin(), models.end(), group.ModelName()) == models.end())
            models.push_back(group.ModelName());
    }

    const bool load_cache = !config.Exists("Cache/Load") || config.GetAs<bool>("Cache/Load");
    const bool save_cache = !config.Exists("Cache/Save") || config.GetAs<bool>("Cache/Save");
    Filesystem::Cache cache;
    if(config.Exists("Cache/Path")) cache.Path() = config.GetAs<std::string>("Cache/Path");
    std::string cache_mode = load_cache && save_cache ? "load and save"
                             : load_cache             ? "load only"
                             : save_cache             ? "save only"
                                                      : "off";
    if(load_cache || save_cache) cache_mode += fmt::format(" ({})", cache.Path().string());

    const std::string seed =
        m_seed_from_card
            ? fmt::format("{} (from the run card)", m_seed)
            : fmt::format("{} (from the clock; set Options/Initialize/Seed to repeat this run)",
                          m_seed);

    std::vector<std::pair<std::string, std::string>> rows;
    rows.emplace_back("Run card", runcard);
    for(size_t i = 0; i < beams.size(); ++i) rows.emplace_back(i == 0 ? "Beam" : "", beams[i]);
    rows.emplace_back("Target", fmt::format("{}", fmt::join(targets, ", ")));
    rows.emplace_back("Nuclear models", fmt::format("{}", fmt::join(models, ", ")));
    rows.emplace_back("Cascade", runCascade ? "on" : "off");
    rows.emplace_back("Hard cuts", doHardCuts ? "on" : "off");
    rows.emplace_back("Events", fmt::format("{}", config["Main/NEvents"].as<size_t>()));
    rows.emplace_back("Output", fmt::format("{} ({})", output, format));
    rows.emplace_back("Seed", seed);
    rows.emplace_back("Cache", cache_mode);
    console::Print(console::Fields("Run configuration", rows));

    if(doHardCuts) PrintCuts();
    if(runCascade) PrintCascade();
}

namespace {

constexpr size_t max_value_width = 90;

// "0.04", "[36.5, 38.5]"
std::string YamlText(const YAML::Node &node) {
    if(!node || node.IsNull()) return "";
    if(node.IsScalar()) return node.as<std::string>();
    YAML::Emitter out;
    out << YAML::Flow << node;
    return out.c_str();
}

// "e- (11), proton (2212)"
std::string ParticleNames(const YAML::Node &pids) {
    std::vector<long> ids;
    if(pids.IsSequence())
        for(const auto &pid : pids) ids.push_back(pid.as<long>());
    else
        ids.push_back(pids.as<long>());

    std::vector<std::string> names;
    for(const auto id : ids) {
        try {
            names.push_back(fmt::format("{} ({})", achilles::ParticleInfo(id).Name(), id));
        } catch(const std::runtime_error &) { names.push_back(fmt::format("{}", id)); }
    }
    return fmt::format("{}", fmt::join(names, ", "));
}

// Nested models by name, lists by length
std::string OptionSummary(const YAML::Node &options) {
    if(!options || options.IsNull()) return "";
    if(!options.IsMap()) return YamlText(options);
    std::vector<std::string> parts;
    for(const auto &option : options) {
        const auto key = option.first.as<std::string>();
        const auto &value = option.second;
        if(value.IsScalar())
            parts.push_back(fmt::format("{}: {}", key, value.as<std::string>()));
        else if(value.IsMap() && value["Name"])
            parts.push_back(fmt::format("{}: {}", key, value["Name"].as<std::string>()));
        else if(value.IsSequence())
            parts.push_back(fmt::format("{}: {} entries", key, value.size()));
        else if(!value.IsNull())
            parts.push_back(key + ": ...");
    }
    return fmt::format("{}", fmt::join(parts, ", "));
}

} // namespace

void achilles::EventGen::PrintCuts() {
    console::Table table({"Cut", "Particles", "Settings"});
    table.SetTitle("Hard cuts");
    for(const auto &cut : config["HardCuts"]) {
        std::vector<std::string> settings;
        for(const auto &setting : cut) {
            const auto key = setting.first.as<std::string>();
            if(key == "Type" || key == "PIDs") continue;
            settings.push_back(fmt::format("{}: {}", key, YamlText(setting.second)));
        }
        table.AddRow(
            {cut["Type"].as<std::string>(), ParticleNames(cut["PIDs"]),
             console::Truncate(fmt::format("{}", fmt::join(settings, ", ")), max_value_width)});
    }
    table.Print();
}

void achilles::EventGen::PrintCascade() {
    console::Table table({"Setting", "Value"});
    table.SetTitle("Cascade");
    for(const auto &setting : config["Cascade"]) {
        const auto key = setting.first.as<std::string>();
        if(key == "Run") continue;
        if(key != "Interactions") {
            table.AddRow({key, console::Truncate(YamlText(setting.second), max_value_width)});
            continue;
        }
        bool first = true;
        for(const auto &interaction : setting.second) {
            std::string value = interaction["Name"] ? interaction["Name"].as<std::string>() : "";
            const auto options = OptionSummary(interaction["Options"]);
            if(!options.empty()) value += " (" + options + ")";
            table.AddRow({first ? "Interactions" : "", console::Truncate(value, max_value_width)});
            first = false;
        }
    }
    table.Print();
}

void achilles::EventGen::PrintProcesses() const {
    using console::Align;
    console::Table table({"Group", "Nucleus", "Nuclear model", "Process"},
                         {Align::right, Align::left, Align::left, Align::left});
    size_t nprocesses = 0;
    for(size_t i = 0; i < process_groups.size(); ++i) {
        const auto &group = process_groups[i];
        if(i != 0) table.AddSeparator();
        bool first = true;
        for(const auto &process : group.Processes()) {
            table.AddRow({first ? group.Id() : "", first ? group.GetNucleus()->ToString() : "",
                          first ? group.ModelName() : "", process.Info().Describe()});
            first = false;
            ++nprocesses;
        }
    }
    table.SetTitle(fmt::format("Processes: {} in {} group{}", nprocesses, process_groups.size(),
                               process_groups.size() == 1 ? "" : "s"));
    table.Print();
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

    // Groups whose results all come from the cache do not need to be saved again
    std::vector<bool> from_cache(process_groups.size(), false);
    if(!config.Exists("Cache/Load") || config.GetAs<bool>("Cache/Load")) {
        for(size_t i = 0; i < process_groups.size(); ++i) {
            auto &group = process_groups[i];
            if(cache.FindCachedState(std::hash<ProcessGroup>{}(group))) {
                spdlog::debug("Loading cached state for {}", group);
                if(!cache.LoadState(group)) {
                    spdlog::warn("Failed to load cached state for {}", group);
                } else {
                    from_cache[i] = group.IsComplete();
                }
            }
        }
    }

    // Generation time is unknown until the unweighting efficiencies are
    auto &tracker = RunTracker::Instance();
    for(const auto &group : process_groups) group.RegisterPendingWork();
    tracker.SetPending("generate", -1);
    tracker.AddStep("generate", "Events");

    spdlog::info("Starting optimization runs for {} process group{}", process_groups.size(),
                 process_groups.size() == 1 ? "" : "s");
    for(size_t i = 0; i < process_groups.size(); ++i) {
        auto &group = process_groups[i];
        group.Optimize();
        m_group_weights.push_back(group.MaxWeight());
        m_max_weight += group.MaxWeight();
        spdlog::debug("Group weights: {} / {}",
                      fmt::join(m_group_weights.begin(), m_group_weights.end(), ", "),
                      m_max_weight);

        if(!from_cache[i] && (!config.Exists("Cache/Save") || config.GetAs<bool>("Cache/Save"))) {
            if(!cache.SaveState(group)) {
                spdlog::warn("Failed to save cached state for {}", group);
            }
        }
    }
    for(auto &wgt : m_group_weights) wgt /= m_max_weight;
    spdlog::debug("Group weights: {}",
                  fmt::join(m_group_weights.begin(), m_group_weights.end(), ", "));

    spdlog::info("Finished optimization");
}

void achilles::EventGen::GenerateEvents(bool batchMode) {
    outputEvents = true;
    console::Configure(batchMode);

    const auto nevents = config["Main/NEvents"].as<size_t>();
    size_t accepted = 0;
    size_t trials = 0;

    auto &tracker = RunTracker::Instance();
    tracker.RemovePending("generate");
    tracker.BeginStage("Event generation");
    tracker.SetStepState("generate", StepState::running);
    spdlog::info("Generating {} events", nevents);

    ProgressBar bar("Generating", static_cast<double>(nevents));
    // Formatting the status is expensive compared to an event
    Stopwatch update_clock;
    auto status = [&]() {
        const double efficiency =
            trials == 0 ? 0 : static_cast<double>(accepted) / static_cast<double>(trials);
        const double elapsed = bar.Elapsed();
        const double rate = elapsed > 0 ? static_cast<double>(accepted) / elapsed : 0;
        const char *sep = console::GetStyle().unicode ? " · " : " | ";
        return fmt::format("{:>5}/{} events{}eff {:4.1f}%{}{:>5}/s",
                           console::FormatCount(static_cast<double>(accepted)),
                           console::FormatCount(static_cast<double>(nevents)), sep,
                           100 * efficiency, sep, console::FormatCount(rate));
    };

    bar.Update(0, status());
    while(accepted < nevents) {
        if(GenerateSingleEvent()) accepted++;
        trials++;
        if(update_clock.Elapsed() >= 0.1) {
            update_clock.Reset();
            bar.Update(static_cast<double>(accepted), status());
        }
    }
    bar.Finish();

    const double efficiency = static_cast<double>(accepted) / static_cast<double>(trials);
    spdlog::info("Generated {} events in {} ({} trials, {:.2f}% efficiency)", accepted,
                 console::FormatDuration(bar.Elapsed()), trials, 100 * efficiency);
    tracker.EndStage(
        fmt::format("{} events, {:.2f}% unweighting efficiency", accepted, 100 * efficiency));
    tracker.SetStepState("generate", StepState::done);

    // Independent integrals: means and variances add
    double xsec = 0, variance = 0;
    for(const auto &group : process_groups) {
        xsec += group.CrossSection().Mean();
        variance += std::pow(group.CrossSection().Error(), 2);
    }
    for(size_t i = 0; i < process_groups.size(); ++i) {
        const auto &group = process_groups[i];
        tracker.AddCrossSection({group.Label(), group.CrossSection().Mean(),
                                 group.CrossSection().Error(), m_group_weights[i],
                                 group.UnweightEff()});
    }
    if(process_groups.size() > 1)
        tracker.AddCrossSection({"Total", xsec, std::sqrt(variance), 1, efficiency});

    const double seconds = bar.Elapsed();
    tracker.SetResult(
        "Events",
        fmt::format("{} written in {} = {} events/day", accepted, console::FormatDuration(seconds),
                    seconds > 0
                        ? console::FormatCount(static_cast<double>(accepted) * 86400 / seconds)
                        : "--"));
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
    // Keep the sampled momentum: this becomes NuHepMC's status-4 incoming beam
    // particle, which consumers read as the probe.
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
    // Cascade the nucleus, with possible re-runs in case of rare errors
    size_t ntrials = 10;
    if(runCascade) {
        spdlog::debug("Runnning cascade");
        for(size_t itrial = 1; itrial <= ntrials; itrial++) {
            try {
                auto tmp_event = event;
                cascade->Evolve(tmp_event, group.GetNucleus());
                event = tmp_event;
                spdlog::trace("Achilles cascade succeeded after {} trials", itrial);
                break;
            } catch(const AchillesCascadeError &e) {
                // Handle rare (~1:1e7) failures from threshold mismatches.
                spdlog::trace("Skipping AchillesCascadeError");
                continue;
            }
            throw AchillesCascadeError("Cascade trials limit reached.");
        }
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
