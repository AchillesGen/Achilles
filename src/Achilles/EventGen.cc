#include "Achilles/EventGen.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Channels.hh"
#include "Achilles/Debug.hh"
#include "Achilles/Event.hh"
#include "Achilles/Exception.hh"
#include "Achilles/Logging.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ReferenceHandler.hh"
#include "Achilles/Units.hh"
#include "Achilles/Writers/EventWriter.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "Plugins/Sherpa/SherpaInterface.hh"
#else
namespace achilles {
class SherpaInterface {};
} // namespace achilles
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdouble-promotion"
#endif
#include "Achilles/Writers/HepMC3EventWriter.hh"
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#include "Achilles/Writers/NuHepMCWriter.hh"

#include "fmt/std.h"
#include "yaml-cpp/yaml.h"

#include <numeric>

achilles::EventGen::EventGen(const std::string &configFile, std::vector<std::string> shargs)
    : config{configFile}, m_run_card{configFile} {
    Setup(std::move(shargs));
}

achilles::EventGen::EventGen(const YAML::Node &node, std::vector<std::string> shargs,
                             const std::string &run_card_name)
    : config{node}, m_run_card{run_card_name} {
    Setup(std::move(shargs));
}

void achilles::EventGen::Setup(std::vector<std::string> shargs) {
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
        if(config.GetAs<int>("Options/Initialize/Seed") > 0)
            seed = config.GetAs<unsigned int>("Options/Initialize/Seed");
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load initial beam and nuclei
    spdlog::trace("Initializing the beams");
    beam = std::make_shared<Beam>(config.GetAs<Beam>("Beam"));
    nuclei = config.GetAs<std::vector<std::shared_ptr<Nucleus>>>("Nuclei");

    // Keep a typed handle only when the configured beam is a geometry beam; the
    // dynamic_cast yields nullptr for every other beam type, leaving the normal
    // (non-geometry) path untouched. All species share one GeometryBeam.
    if(!beam->BeamIDs().empty())
        m_geometry_beam = dynamic_cast<GeometryBeam *>(beam->at(*beam->BeamIDs().begin()).get());

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
        spdlog::info("Apply hard cuts? {}", doHardCuts);
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
    } else if(format == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(name, zipped);
    } else if(format == "NuHepMC") {
        writer = std::make_unique<NuHepMCWriter>(name, zipped);
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}", format);
        throw std::runtime_error(msg);
    }
    writer->WriteHeader(m_run_card, process_groups);
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
                spdlog::debug("Loading cached state for {}", group);
                if(!cache.LoadState(group)) {
                    spdlog::warn("Failed to load cached state for {}", group);
                }
            }
        }
    }

    spdlog::info("Starting optimization runs");
    for(auto &group : process_groups) {
        group.SetupSplines();
        group.Optimize();
        m_group_weights.push_back(group.MaxWeight());
        m_max_weight += group.MaxWeight();
        spdlog::info("Group weights: {} / {}",
                     fmt::join(m_group_weights.begin(), m_group_weights.end(), ", "), m_max_weight);
        std::cout << "Estimated unweighting eff for this group: ";
        for(auto &process : group.Processes()) {
            std::cout << process.UnweightEff() << " ";
            std::cout << std::endl;
        }

        if(!config.Exists("Cache/Save") || config.GetAs<bool>("Cache/Save")) {
            if(!cache.SaveState(group)) {
                spdlog::warn("Failed to save cached state for {}", group);
            }
        }
    }
    for(auto &wgt : m_group_weights) wgt /= m_max_weight;
    spdlog::info("Group weights: {}",
                 fmt::join(m_group_weights.begin(), m_group_weights.end(), ", "));
    spdlog::info("Finished optimization");
}

void achilles::EventGen::GenerateEvents(bool batchMode) {
    outputEvents = true;

    const auto nevents = config["Main/NEvents"].as<size_t>();
    size_t accepted = 0;
    size_t statusUpdate = 1;
    size_t lastUpdate = 0; // Prevents the same # of events from being logged more than once
                           // (would happen when events were rejected)

    auto spdlog_info = [](size_t acc, size_t nEv) {
        spdlog::info("Generated {} / {} events", acc, nEv);
    };
    auto fmt_print = [](size_t acc, size_t nEv) {
        fmt::print("Generated {} / {} events\r", acc, nEv);
    };
    auto printFormat = batchMode ? spdlog_info : fmt_print;

    printFormat(0, nevents);
    while(accepted < nevents) {
        if(accepted % statusUpdate == 0 && accepted > lastUpdate) {
            printFormat(accepted, nevents);
            lastUpdate = accepted;
            if(accepted >= 10 * statusUpdate) statusUpdate *= 10;
        }
        if(GenerateSingleEvent()) accepted++;
    }
    printFormat(accepted, nevents);
}

void achilles::EventGen::InjectRay(const FourVector &lab_ray, PID nu, PID target) {
    if(!m_geometry_beam)
        throw std::runtime_error("EventGen::InjectRay requires a beam of Type: Geometry");
    m_geometry_beam->SetMode(GeometryBeam::Mode::Generate);
    m_geometry_beam->SetInjected(lab_ray);
    m_injection = std::make_pair(nu, target);
    spdlog::trace("EventGen::InjectRay: E = {} MeV, nu = {}, target = {}", lab_ray.E(), nu, target);
}

void achilles::EventGen::WriteRunInfo(std::shared_ptr<HepMC3::GenRunInfo> run) const {
    NuHepMCWriter::WriteHeader(std::move(run), m_run_card, process_groups);
}

void achilles::EventGen::Finalize() {
    // TODO: flush splines / run statistics as the geometry workflow requires.
    spdlog::info("EventGen: finalizing run");
}

double achilles::EventGen::VertexEnvelope(PID nu, PID target) const {
    const auto weights = GeometryGroupWeights(process_groups, nu, target);
    // Internal weights are in nb; NuGeom works in cm^2.
    const double envelope =
        std::accumulate(weights.begin(), weights.end(), 0.0) * Constant::NB_TO_CM2;
    spdlog::trace("EventGen::VertexEnvelope(nu = {}, target = {}) = {} cm^2", nu, target, envelope);
    return envelope;
}

double achilles::EventGen::TotalXSec(double energy, PID nu, PID target) const {
    double total = 0;
    for(const auto &group : process_groups)
        if(group.GetNucleus()->ID() == target) total += group.TotalXSec(energy, nu);
    // Spline values are in nb (same units as the raw weights); convert to cm^2
    // so both layer-1 accessors feed NuGeom in the same units.
    total *= Constant::NB_TO_CM2;
    spdlog::trace("EventGen::TotalXSec(E = {} MeV, nu = {}, target = {}) = {} cm^2", energy, nu,
                  target, total);
    return total;
}

bool achilles::EventGen::GenerateSingleEvent() {
    // Select the process group and generate an event. In geometry mode the
    // injected ray restricts selection to the matching (neutrino, target);
    // groups are weighted by only the matching species' max-weight so a group's
    // unrelated (other-species) processes don't bias the choice. The normal path
    // (no injection) is unchanged.
    std::optional<PID> required_nu;
    size_t group_idx;
    if(m_injection) {
        const auto nu = m_injection->first;
        const auto target = m_injection->second;
        required_nu = nu;

        // The warm-up flat flux only calibrated the integrator and unweighter
        // inside the configured EnergyRange; rays outside it have no model
        // support and must be rejected (zero cross section), not extrapolated.
        const double e_inj = m_geometry_beam->Injected().E();
        if(e_inj < m_geometry_beam->MinEnergy() || e_inj > m_geometry_beam->MaxEnergy()) {
            if(m_range_warnings < kMaxRangeWarnings) {
                spdlog::warn("EventGen: injected ray energy {} MeV outside beam EnergyRange "
                             "[{}, {}]; treating as zero cross section. Widen EnergyRange to "
                             "cover the flux file. ({} more warnings before silencing)",
                             e_inj, m_geometry_beam->MinEnergy(), m_geometry_beam->MaxEnergy(),
                             kMaxRangeWarnings - m_range_warnings - 1);
                ++m_range_warnings;
            } else {
                spdlog::debug("EventGen: injected ray energy {} MeV outside beam EnergyRange; "
                              "rejected",
                              e_inj);
            }
            return false;
        }

        auto weights = GeometryGroupWeights(process_groups, nu, target);
        if(std::accumulate(weights.begin(), weights.end(), 0.0) == 0.0) {
            spdlog::debug("EventGen: no process available for nu = {} on target = {}; rejected", nu,
                          target);
            return false;
        }
        group_idx = Random::Instance().SelectIndex(weights);
        spdlog::trace("EventGen: selected process group {} for nu = {}, target = {}", group_idx, nu,
                      target);
    } else {
        group_idx = Random::Instance().SelectIndex(m_group_weights);
    }
    auto &group = process_groups[group_idx];
    auto &&event = group.GenerateEvent(required_nu);
    // Weight convention depends on the run mode:
    //  - Analytic (no injected ray): rescale the unweighted (unit) event back to
    //    sigma units so sum(weight)/N_trials = sigma_tot. m_max_weight is the
    //    global all-group envelope used by the analytic estimator.
    //  - Geometry / EnvelopeNoRetry: leave the partial-unweight result
    //    max(1, raw_w/w_max_p); the max-weight scale is absorbed into NuGeom's
    //    layer-1 envelope (= VertexEnvelope), and the rate normalization is
    //    carried by NuGeom's POT / 1-over-max_prob bookkeeping.
    //  - Geometry / TotalXSecRetry: NuGeom's layer-1 uses sigma_tot, so the
    //    event carries the LHC convention max(w_max_p, raw_w) = w_max_p * Unweight.
    // Geometry mode (m_injection): the event keeps the partial-unweight
    // weight max(1, raw_w/w_max_p) ~ 1 in BOTH sampling modes (E.C.1 "CV"
    // convention); the sigma scale for the driver's estimator is
    // VertexEnvelope, applied driver-side per trial.
    if(!m_injection) event.Weight() *= m_max_weight;
    if(event.Weight() == 0) {
        spdlog::trace("EventGen: trial rejected by unweighter (weight = 0)");
        // In geometry mode NuGeom owns the output and treats false as a rejected
        // ray, so don't write anything here.
        if(!m_injection) writer->Write(event);
        return false;
    }
    spdlog::trace("EventGen: trial accepted with weight = {}", event.Weight());

    // Geometry mode: undo the rotation that aligned the injected ray with +z
    // before anything else, so the whole event -- history and cascade included --
    // is built and stored in the lab frame.
    if(m_injection) {
        const auto &rotation = m_geometry_beam->Rotation();
        for(auto &lepton : event.Leptons()) rotation.RotateBack(lepton.Momentum());
        for(auto &hadron : event.Hadrons()) rotation.RotateBack(hadron.Momentum());
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
    // Setup beam in history. The beam particle is the incoming neutrino itself
    // (leptons[0], the initial-state lepton at its sampled energy/direction), not
    // a hardcoded max-energy +z vector.
    auto init_lep = event.Leptons()[0];
    auto init_beam = init_lep;
    init_beam.Status() = ParticleStatus::beam;
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

    // In geometry mode, hand the event back to NuGeom (via LastEvent) instead of
    // writing it ourselves; it is already in the lab frame.
    if(m_injection) {
        m_last_event = event;
        return true;
    }

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
