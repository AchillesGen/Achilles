#include "Achilles/Process.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Channels.hh"
#include "Achilles/Exceptions.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/MultiChannel.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Settings.hh"
#include "Achilles/Unweighter.hh"
#include "Achilles/Utilities.hh"

#include "fmt/ranges.h"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/SherpaInterface.hh"
#else
// Create dummy class
namespace achilles {
class SherpaInterface {};
} // namespace achilles
#endif

using achilles::Process;
using achilles::ProcessGroup;
using achilles::refParticles;

refParticles Process::SelectParticles(const refParticles &protons, const refParticles &neutrons,
                                      const std::vector<PID> &ids,
                                      const std::vector<FourVector> &momenta,
                                      ParticleStatus status) const {
    refParticles result;
    for(size_t i = 0; i < ids.size(); ++i) {
        if(ids[i] == PID::proton()) {
            result.push_back(protons[i]);
        } else {
            result.push_back(neutrons[i]);
        }

        auto &part = result.back().get();
        part.Momentum() = momenta[i];
        part.Status() = status;
    }

    return result;
}

void Process::SetupHadrons(Event &event) const {
    FourVector lep_in;
    std::vector<FourVector> had_in, lep_out, had_out, spect;
    ExtractMomentum(event, lep_in, had_in, lep_out, had_out, spect);
    std::vector<Particle> leptons, hadrons;

    // Setup hadrons
    if(had_in.size() != m_info.m_hadronic.first.size() ||
       had_out.size() != m_info.m_hadronic.second.size())
        throw std::runtime_error(
            "Process: Number of hadrons momenta does not match number of hadrons in event");

    // Handle coherent scattering as a special case
    if(ParticleInfo(m_info.m_hadronic.first[0]).IsNucleus()) {
        event.Hadrons().clear();
        event.Hadrons().resize(2);
        const auto id_in = m_info.m_hadronic.first[0];
        const auto id_out = m_info.m_hadronic.second[0];

        Particle initial(id_in, had_in[0]);
        initial.Status() = ParticleStatus::initial_state;
        event.Hadrons()[0] = initial;

        Particle final(id_out, had_out[0]);
        final.Status() = ParticleStatus::final_state;
        event.Hadrons()[1] = final;

        return;
    }

    // TODO: Handle hydrogen

    // NOTE: The initial state and spectator have to be handled separately
    //       to ensure that they aren't the same particle

    // Select initial state nucleons
    auto protons = event.Protons(ParticleStatus::background);
    auto neutrons = event.Neutrons(ParticleStatus::background);
    size_t nsamples = m_info.m_hadronic.first.size();
    auto sampled_protons = Random::Instance().Sample(nsamples, protons);
    auto sampled_neutrons = Random::Instance().Sample(nsamples, neutrons);

    auto initial_states =
        SelectParticles(sampled_protons, sampled_neutrons, m_info.m_hadronic.first, had_in,
                        ParticleStatus::initial_state);
    spdlog::debug("Selected initial states: {}", fmt::join(initial_states, ", "));

    // Sample for spectators
    protons = event.Protons(ParticleStatus::background);
    neutrons = event.Neutrons(ParticleStatus::background);
    nsamples = m_info.m_hadronic.first.size();
    sampled_protons = Random::Instance().Sample(nsamples, protons);
    sampled_neutrons = Random::Instance().Sample(nsamples, neutrons);
    auto spectators = SelectParticles(sampled_protons, sampled_neutrons, m_info.m_spectator, spect,
                                      ParticleStatus::spectator);
    spdlog::debug("Selected spectators states: {}", fmt::join(initial_states, ", "));

    // Initialize final state hadrons
    // TODO: Handle propagating deltas
    // TODO: Handle selecting position for things like MEC+pion production
    size_t cur_idx = 0;
    ThreeVector position;
    for(size_t i = 0; i < had_out.size(); ++i) {
        Particle part(m_info.m_hadronic.second[i]);
        if(ParticleInfo(m_info.m_hadronic.second[i]).IsBaryon()) {
            position = initial_states[cur_idx++].get().Position();
        }
        part.Status() = ParticleStatus::propagating;
        part.Momentum() = had_out[i];
        part.Position() = position;
        event.Hadrons().push_back(part);
    }
}

void Process::ExtractMomentum(const Event &event, FourVector &lep_in,
                              std::vector<FourVector> &had_in, std::vector<FourVector> &lep_out,
                              std::vector<FourVector> &had_out,
                              std::vector<FourVector> &spect) const {
    const auto &momentum = event.Momentum();
    if(m_info.m_mom_map.size() == 0) {
        static constexpr size_t lepton_in_end_idx = 1;
        const size_t hadron_end_idx = m_info.m_hadronic.first.size() + lepton_in_end_idx;
        const size_t lepton_end_idx = m_info.m_leptonic.second.size() + hadron_end_idx;
        const size_t had_out_end_idx = m_info.m_hadronic.second.size() + lepton_end_idx;

        spdlog::trace("{}, {}, {}, {}, {}", lepton_in_end_idx, hadron_end_idx, lepton_end_idx,
                      had_out_end_idx, momentum.size());
        spdlog::trace("{}", m_info);

        lep_in = momentum[0];
        had_in = std::vector<FourVector>(momentum.begin() + lepton_in_end_idx,
                                         momentum.begin() + static_cast<int>(hadron_end_idx));
        lep_out = std::vector<FourVector>(momentum.begin() + static_cast<int>(hadron_end_idx),
                                          momentum.begin() + static_cast<int>(lepton_end_idx));
        had_out = std::vector<FourVector>(momentum.begin() + static_cast<int>(lepton_end_idx),
                                          momentum.begin() + static_cast<int>(had_out_end_idx));
        spect = std::vector<FourVector>(momentum.begin() + static_cast<int>(had_out_end_idx),
                                        momentum.end());
    } else {
        lep_in = momentum[0];
        had_in.push_back(momentum[1]);
        for(size_t i = 2; i < m_info.m_mom_map.size(); ++i) {
            PID pid(m_info.m_mom_map[i]);
            bool is_lepton = ParticleInfo(pid).IsLepton();
            if(is_lepton)
                lep_out.push_back(momentum[i]);
            else
                had_out.push_back(momentum[i]);
        }
    }
}

void Process::ExtractParticles(const Event &event, Particle &lep_in, std::vector<Particle> &had_in,
                               std::vector<Particle> &lep_out, std::vector<Particle> &had_out,
                               std::vector<Particle> &spect) const {
    spdlog::trace("Extracting particles from event");
    const auto &momentum = event.Momentum();
    if(m_info.m_mom_map.size() == 0) {
        static constexpr size_t lepton_in_end_idx = 1;
        const size_t hadron_end_idx = m_info.m_hadronic.first.size() + lepton_in_end_idx;
        const size_t lepton_end_idx = m_info.m_leptonic.second.size() + hadron_end_idx;
        const size_t had_out_end_idx = m_info.m_hadronic.second.size() + lepton_end_idx;
        lep_in = Particle(m_info.m_leptonic.first, momentum[0]);
        for(size_t i = 0; i < m_info.m_hadronic.first.size(); ++i) {
            had_in.emplace_back(m_info.m_hadronic.first[i], momentum[i + lepton_in_end_idx]);
        }
        for(size_t i = 0; i < m_info.m_leptonic.second.size(); ++i) {
            lep_out.emplace_back(m_info.m_leptonic.second[i], momentum[i + hadron_end_idx]);
        }
        for(size_t i = 0; i < m_info.m_hadronic.second.size(); ++i) {
            had_out.emplace_back(m_info.m_hadronic.second[i], momentum[i + lepton_end_idx]);
        }
        for(size_t i = 0; i < m_info.m_spectator.size(); ++i) {
            spect.emplace_back(m_info.m_spectator[i], momentum[i + had_out_end_idx]);
        }
    } else {
        lep_in = Particle(m_info.m_mom_map[1], momentum[0]);
        had_in.emplace_back(m_info.m_mom_map[0], momentum[1]);
        for(size_t i = 2; i < m_info.m_mom_map.size(); ++i) {
            PID pid(m_info.m_mom_map[i]);
            bool is_lepton = ParticleInfo(pid).IsLepton();
            spdlog::trace("Extracting particle {}, with PID: {}", i, pid);
            if(is_lepton)
                lep_out.emplace_back(pid, momentum[i]);
            else
                had_out.emplace_back(pid, momentum[i]);
        }
    }
}

achilles::FourVector Process::ExtractQ(const Event &event) const {
    const auto &momentum = event.Momentum();
    auto q = momentum[0];
    size_t nleptons = m_info.m_leptonic.second.size();
    for(size_t i = 0; i < nleptons; ++i) { q -= momentum[i + 1 + m_info.m_hadronic.first.size()]; }
    return q;
}

void Process::SetID(achilles::NuclearModel *model) {
    // Get range for interaction type
    m_id = ToID(model->Mode());

    // Get range for CC or NC
    auto charge = m_info.LeptonicCharge();
    if(std::abs(charge) == 0) { m_id += 50; }

    // Add number of protons in initial state
    for(const auto &nucleons : m_info.m_hadronic.first)
        if(nucleons == PID::proton()) m_id++;

    // Add number of pions in the final state
    for(const auto &particle : m_info.m_hadronic.second)
        if(particle == PID::pion0() || std::abs(particle.AsInt()) == PID::pionp().AsInt()) m_id++;
}

std::string Process::Name(achilles::XSecBackend *backend) const {
    // TODO: Clean this up to better handle the header writing for cascade mode
    if(m_info.m_leptonic.first == PID::undefined()) return "Cascade";

    std::string name = backend->GetNuclearModel()->GetName();
    name += (m_id % 100) / 50 == 0 ? "CC" : "NC";
    // Add number of protons in initial state
    int nprotons = 0;
    for(const auto &nucleons : m_info.m_hadronic.first)
        if(nucleons == PID::proton()) nprotons++;
    // Add number of pions in the final state
    int npions = 0;
    for(const auto &particle : m_info.m_hadronic.second)
        if(particle == PID::pion0() || std::abs(particle.AsInt()) == PID::pionp().AsInt()) npions++;
    name += std::to_string(nprotons) + "p";
    name += std::to_string(npions) + "pi";
    return name;
}

std::string Process::Description(achilles::XSecBackend *backend) const {
    // TODO: Clean this up to better handle the header writing for cascade mode
    if(m_info.m_leptonic.first == PID::undefined()) return "Cascade Only";

    std::string description = backend->GetNuclearModel()->GetName() + " ";
    std::stringstream ss;
    ss << m_info;
    description += ss.str();

    return description;
}

std::string Process::InspireHEP(achilles::XSecBackend *backend) const {
    // TODO: Clean this up to better handle the header writing for cascade mode
    if(m_info.m_leptonic.first == PID::undefined()) return "";

    return backend->GetNuclearModel()->InspireHEP();
}

achilles::ProcessMetadata Process::Metadata(achilles::XSecBackend *backend) const {
    return {ID(), Name(backend), Description(backend), InspireHEP(backend)};
}

bool Process::SaveState(std::ostream &os) const {
    m_info.SaveState(os);
    m_xsec.SaveState(os);
    m_unweighter->SaveState(os);
    return true;
}

bool Process::LoadState(std::istream &is) {
    m_info.LoadState(is);
    m_xsec.LoadState(is);
    m_unweighter->LoadState(is);
    return true;
}

std::vector<int> ProcessGroup::ProcessIds() const {
    std::vector<int> data;
    for(auto &process : m_processes) { data.push_back(process.ID()); }
    return data;
}

std::vector<achilles::ProcessMetadata> ProcessGroup::Metadata() const {
    std::vector<achilles::ProcessMetadata> data;
    for(auto &process : m_processes) { data.push_back(process.Metadata(m_backend.get())); }
    return data;
}

void ProcessGroup::SetupLeptons(Event &event, std::optional<size_t> process_idx) const {
    spdlog::trace("Setting up leptons");
    FourVector lep_in;
    std::vector<FourVector> had_in, lep_out, had_out, spect;
    auto &process = m_processes[process_idx.value_or(0)];
    process.ExtractMomentum(event, lep_in, had_in, lep_out, had_out, spect);
    std::vector<Particle> leptons;

    // Setup leptons
    const auto &info = process.Info();
    leptons.emplace_back(info.m_leptonic.first, lep_in);
    leptons.back().Status() = ParticleStatus::initial_state;
    if(lep_out.size() != info.m_leptonic.second.size())
        throw std::runtime_error(
            "Process: Number of lepton momenta does not match number of leptons in event");
    for(size_t i = 0; i < lep_out.size(); ++i) {
        leptons.emplace_back(info.m_leptonic.second[i], lep_out[i]);
        leptons.back().Status() = ParticleStatus::final_state;
    }

    event.Leptons() = leptons;
}

void ProcessGroup::CrossSection(Event &event, std::optional<size_t> process_idx) {
    if(!process_idx) {
        double weight = 0;
        for(size_t i = 0; i < m_processes.size(); ++i) {
            auto process_weight = m_backend->CrossSection(event, m_processes[i]);
            if(b_calc_weights) m_processes[i].AddWeight(process_weight);
            weight += process_weight;
        }
        event.Weight() = weight;
        if(b_calc_weights) m_xsec += weight;
    } else {
        auto &process = m_processes[process_idx.value()];
        auto weight = m_backend->CrossSection(event, process);
        event.Weight() = process.Unweight(weight);
    }
}

size_t ProcessGroup::SelectProcess() const {
    return Random::Instance().SelectIndex(m_process_weights);
}

std::map<size_t, ProcessGroup> ProcessGroup::ConstructGroups(const Settings &settings,
                                                             NuclearModel *model,
                                                             std::shared_ptr<Beam> beam,
                                                             std::shared_ptr<Nucleus> nucleus) {
    std::map<size_t, ProcessGroup> groups;
    for(const auto &process_node : settings["Processes"]) {
        auto process_info = process_node.as<ProcessInfo>();
        auto infos = model->AllowedStates(process_info);
        const auto unweight_name = settings.GetAs<std::string>("Options/Unweighting/Name");
        for(auto &info : infos) {
            try {
                auto unweighter =
                    UnweighterFactory::Initialize(unweight_name, settings["Options/Unweighting"]);
                Process process(info, std::move(unweighter));
                process.SetID(model);
                const auto multiplicity = info.Multiplicity();
                if(groups.find(multiplicity) == groups.end()) {
                    groups.insert({multiplicity, ProcessGroup(beam, nucleus)});
                }
                groups.at(multiplicity).AddProcess(std::move(process));
            } catch(std::out_of_range &e) {
                spdlog::error("Unweighter: Requested unweighter \"{}\", did you mean \"{}\"",
                              unweight_name,
                              achilles::GetSuggestion(UnweighterFactory::List(), unweight_name));
                exit(-1);
            }
        }
    }

    return groups;
}

void ProcessGroup::SetupBackend(const Settings &settings, std::unique_ptr<NuclearModel> model,
                                SherpaInterface *sherpa) {
    auto backend_name = settings.GetAs<std::string>("Backend/Name");
    m_backend = XSecBuilder(backend_name)
                    .AddOptions(settings["Backend/Options"])
                    .AddNuclearModel(std::move(model))
                    .AddSherpa(sherpa)
                    .build();

    for(auto &process : m_processes) m_backend->AddProcess(process);
}

bool ProcessGroup::SetupIntegration(const Settings &config) {
    if(m_processes.size() == 0)
        throw std::runtime_error("ProcessGroup: Number of processes found is zero!");

    try {
        m_backend->SetupChannels(m_processes[0].Info(), m_beam, m_integrand, m_nucleus->ID());
    } catch(const InvalidChannel &) { return false; }

    MultiChannelParams multichannel_params;
    if(config.Exists("Options/Initialize/Parameters"))
        multichannel_params = config.GetAs<MultiChannelParams>("Options/Initialize/Parameters");
    m_integrator = MultiChannel(m_integrand.NDims(), m_integrand.NChannels(), multichannel_params);
    if(config.Exists("Options/Initialize/Accuracy"))
        m_integrator.Parameters().rtol = config.GetAs<double>("Options/Initialize/Accuracy");

    return true;
}

bool ProcessGroup::NeedsOptimization() const {
    if(!m_integrator.HasResults()) return true;
    double rel_err = m_integrator.LastResult().RelError();
    return m_integrator.NeedsOptimization(rel_err);
}

void ProcessGroup::Optimize() {
    // TODO: Clean up this control flow. It is a bit of a mess
    b_optimize = true;

    auto func = [&](const std::vector<FourVector> &mom, const double &wgt) {
        return SingleEvent(mom, wgt).Weight();
    };
    m_integrand.Function() = func;
    if(NeedsOptimization()) {
        spdlog::info(
            "Optimizing process group: Nucleus = {}, Nuclear Model = {}, Multiplicity = {}",
            m_nucleus->ToString(), m_backend->GetNuclearModel()->GetName(),
            m_processes[0].Info().Multiplicity());
        m_integrator.Optimize(m_integrand);
    }

    // Ensure that the integrator does not save these results into the summary
    m_integrator.UpdateSummary(false);
    if(m_maxweight == 0) {
        spdlog::info("Calculating maximum weight");
        b_calc_weights = true;
        m_integrator.Parameters().ncalls = 100000;
        for(size_t i = 0; i < 3; ++i) m_integrator(m_integrand);
        b_calc_weights = false;

        // Store max weight and weight vector
        m_process_weights.resize(m_processes.size());
        for(size_t i = 0; i < m_processes.size(); ++i) {
            m_process_weights[i] = m_processes[i].MaxWeight();
            m_maxweight += m_process_weights[i];
            spdlog::info("Process xsec: {} ", m_processes[i].TotalCrossSection());
        }
    } else {
        for(size_t i = 0; i < m_processes.size(); ++i) {
            spdlog::info("Process xsec: {} ", m_processes[i].TotalCrossSection());
        }
    }
    b_optimize = false;
    spdlog::info("Total xsec: {} +/- {} ({}%)", m_xsec.Mean(), m_xsec.Error(),
                 m_xsec.Error() / m_xsec.Mean() * 100);
    spdlog::info("Process weights: {} / {}",
                 fmt::join(m_process_weights.begin(), m_process_weights.end(), ", "), m_maxweight);
    for(size_t i = 0; i < m_processes.size(); ++i) { m_process_weights[i] /= m_maxweight; }
}

void ProcessGroup::Summary() const {}

achilles::Event ProcessGroup::GenerateEvent() {
    const auto mom = m_integrator.GeneratePoint(m_integrand);
    const auto ps_wgt = m_integrator.GenerateWeight(m_integrand, mom);
    return SingleEvent(mom, ps_wgt);
}

achilles::Event ProcessGroup::SingleEvent(const std::vector<FourVector> &mom, double ps_wgt) {
    Event event = Event(m_nucleus, mom, ps_wgt);

    spdlog::debug("Event Phase Space:");
    // size_t idx = 0;
    // for(const auto &momentum : event.Momentum()) {
    //     spdlog::debug("\t{}: {} (M2 = {})", ++idx, momentum, momentum.M2());
    // }
    // Cut on leptons: NOTE: This assumes that all processes in the group have the same leptons
    auto process_opt = b_optimize ? std::nullopt : std::optional<size_t>(SelectProcess());
    SetupLeptons(event, process_opt);
    if(!m_cuts.EvaluateCuts(event.Particles())) {
        // Ensure process weights are tracked correctly
        if(b_calc_weights) {
            for(auto &process : m_processes) process.AddWeight(0);
            m_xsec += 0;
        }
        event.Weight() = 0;
        return event;
    }

    CrossSection(event, process_opt);

    // If training the integrator or weight is zero, we can stop here
    if(b_optimize || event.Weight() == 0) return event;

    // Otherwise, we need to fill the event with the selected process
    auto &process = m_processes[process_opt.value()];
    event.Flux() = m_beam->EvaluateFlux(process.Info().m_leptonic.first, event.Momentum()[0]);
    event.ProcessId() = process.ID();
    process.SetupHadrons(event);

    return event;
}

std::vector<int> achilles::AllProcessIDs(const std::vector<ProcessGroup> &groups) {
    std::vector<int> ids;
    for(const auto &group : groups) {
        auto group_ids = group.ProcessIds();
        ids.insert(ids.end(), group_ids.begin(), group_ids.end());
    }

    return ids;
}

std::vector<achilles::ProcessMetadata>
achilles::AllProcessMetadata(const std::vector<ProcessGroup> &groups) {
    std::vector<ProcessMetadata> data;
    for(const auto &group : groups) {
        auto group_data = group.Metadata();
        data.insert(data.end(), group_data.begin(), group_data.end());
    }

    return data;
}

bool ProcessGroup::Save(const fs::path &cache_dir) const {
    std::ofstream out_processes(cache_dir / "processes.achilles");
    out_processes << m_processes.size() << " ";
    for(const auto &process : m_processes) process.SaveState(out_processes);

    std::ofstream out_integrator(cache_dir / "integrator.achilles");
    m_integrator.SaveState(out_integrator);
    m_integrand.SaveState(out_integrator);

    std::ofstream out_xsec(cache_dir / "xsec.achilles");
    m_xsec.SaveState(out_xsec);

    // TODO: Store some metadata for validation
    std::ofstream out_metadata(cache_dir / "metadata.achilles");

    std::ofstream out_state(cache_dir / "state.achilles");
    out_state << b_optimize << " " << b_calc_weights << " " << m_maxweight << " ";
    out_state << m_process_weights.size() << " ";
    for(const auto &weight : m_process_weights) out_state << weight << " ";

    return true;
}

// TODO: Add validation checks
bool ProcessGroup::Load(const fs::path &cache_dir) {
    std::ifstream in_processes(cache_dir / "processes.achilles");
    size_t nprocesses;
    in_processes >> nprocesses;
    if(nprocesses != m_processes.size())
        throw std::runtime_error("ProcessGroup: Number of processes does not match cache");
    for(auto &process : m_processes) process.LoadState(in_processes);

    std::ifstream in_integrator(cache_dir / "integrator.achilles");
    m_integrator.LoadState(in_integrator);
    m_integrand.LoadState(in_integrator);

    std::ifstream in_xsec(cache_dir / "xsec.achilles");
    m_xsec.LoadState(in_xsec);

    // TODO: Load some metadata for validation
    std::ifstream in_metadata(cache_dir / "metadata.achilles");

    std::ifstream in_state(cache_dir / "state.achilles");
    in_state >> b_optimize >> b_calc_weights >> m_maxweight;
    size_t nweights;
    in_state >> nweights;
    m_process_weights.resize(nweights);
    for(auto &weight : m_process_weights) in_state >> weight;

    return true;
}

std::size_t std::hash<ProcessGroup>::operator()(const achilles::ProcessGroup &p) const {
    std::size_t seed = 0;
    for(const auto &process : p.Processes()) {
        std::hash<achilles::Process> hasher;
        seed ^= hasher(process) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    seed ^= std::hash<achilles::Beam>{}(*(p.m_beam.get())) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<achilles::Nucleus>{}(*(p.m_nucleus.get())) + 0x9e3779b9 + (seed << 6) +
            (seed >> 2);
    return seed;
}
