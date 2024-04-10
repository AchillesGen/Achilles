#include "Achilles/Process.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Channels.hh"
#include "Achilles/Exceptions.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

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

    // Select initial state nucleons and spectators
    auto protons = event.Protons(ParticleStatus::background);
    auto neutrons = event.Neutrons(ParticleStatus::background);
    size_t nsamples = m_info.m_hadronic.first.size() + m_info.m_spectator.size();
    auto sampled_protons = Random::Instance().Sample(nsamples, protons);
    auto sampled_neutrons = Random::Instance().Sample(nsamples, neutrons);

    auto initial_states =
        SelectParticles(sampled_protons, sampled_neutrons, m_info.m_hadronic.first, had_in,
                        ParticleStatus::initial_state);
    auto spectators = SelectParticles(sampled_protons, sampled_neutrons, m_info.m_spectator, spect,
                                      ParticleStatus::spectator);

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
    static constexpr size_t lepton_in_end_idx = 1;
    const size_t hadron_end_idx = m_info.m_hadronic.first.size() + lepton_in_end_idx;
    const size_t lepton_end_idx = m_info.m_leptonic.second.size() + hadron_end_idx;
    const size_t had_out_end_idx = m_info.m_hadronic.second.size() + lepton_end_idx;
    const auto &momentum = event.Momentum();
    lep_in = momentum[0];
    had_in = std::vector<FourVector>(momentum.begin() + lepton_in_end_idx,
                                     momentum.begin() + static_cast<int>(hadron_end_idx));
    lep_out = std::vector<FourVector>(momentum.begin() + static_cast<int>(hadron_end_idx),
                                      momentum.begin() + static_cast<int>(lepton_end_idx));
    had_out = std::vector<FourVector>(momentum.begin() + static_cast<int>(lepton_end_idx),
                                      momentum.begin() + static_cast<int>(had_out_end_idx));
    spect = std::vector<FourVector>(momentum.begin() + static_cast<int>(had_out_end_idx),
                                    momentum.end());
}

void Process::ExtractParticles(const Event &event, Particle &lep_in, std::vector<Particle> &had_in,
                               std::vector<Particle> &lep_out, std::vector<Particle> &had_out,
                               std::vector<Particle> &spect) const {
    static constexpr size_t lepton_in_end_idx = 1;
    const size_t hadron_end_idx = m_info.m_hadronic.first.size() + lepton_in_end_idx;
    const size_t lepton_end_idx = m_info.m_leptonic.second.size() + hadron_end_idx;
    const size_t had_out_end_idx = m_info.m_hadronic.second.size() + lepton_end_idx;
    const auto &momentum = event.Momentum();
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
    for(size_t i = 0; i < spect.size(); ++i) {
        spect.emplace_back(m_info.m_spectator[i], momentum[i + had_out_end_idx]);
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
    std::string description = backend->GetNuclearModel()->GetName() + " ";
    std::stringstream ss;
    ss << m_info;
    description += ss.str();

    return description;
}

std::string Process::InspireHEP(achilles::XSecBackend *backend) const {
    return backend->GetNuclearModel()->InspireHEP();
}

achilles::ProcessMetadata Process::Metadata(achilles::XSecBackend *backend) const {
    return {ID(), Name(backend), Description(backend), ""};
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

std::map<size_t, ProcessGroup> ProcessGroup::ConstructGroups(const YAML::Node &node,
                                                             NuclearModel *model,
                                                             std::shared_ptr<Beam> beam,
                                                             std::shared_ptr<Nucleus> nucleus) {
    std::map<size_t, ProcessGroup> groups;
    for(const auto &process_node : node["Processes"]) {
        auto process_info = process_node.as<ProcessInfo>();
        auto infos = model->AllowedStates(process_info);
        const auto unweight_name = node["Unweighting"]["Name"].as<std::string>();
        for(auto &info : infos) {
            auto unweighter = UnweighterFactory::Initialize(unweight_name, node["Unweighting"]);
            Process process(info, std::move(unweighter));
            process.SetID(model);
            const auto multiplicity = info.Multiplicity();
            if(groups.find(multiplicity) == groups.end()) {
                groups.insert({multiplicity, ProcessGroup(beam, nucleus)});
            }
            groups.at(multiplicity).AddProcess(std::move(process));
        }
    }

    return groups;
}

void ProcessGroup::SetupBackend(const YAML::Node &node, std::unique_ptr<NuclearModel> model,
                                SherpaInterface *sherpa) {
    auto backend_name = node["Backend"]["Name"].as<std::string>();
    m_backend = XSecBuilder(backend_name)
                    .AddOptions(node["Backend"]["Options"])
                    .AddNuclearModel(std::move(model))
                    .AddSherpa(sherpa)
                    .build();

    for(auto &process : m_processes) m_backend->AddProcess(process);
}

bool ProcessGroup::SetupIntegration(const YAML::Node &config) {
    if(m_processes.size() == 0)
        throw std::runtime_error("ProcessGroup: Number of processes found is zero!");

    try {
        m_backend->SetupChannels(m_processes[0].Info(), m_beam, m_integrand, m_nucleus->ID());
    } catch(const InvalidChannel &) { return false; }

    // TODO: Fix scaling to be consistent with Chili paper
    m_integrator = MultiChannel(m_integrand.NDims(), m_integrand.NChannels(), {1000, 2});
    if(config["Initialize"]["Accuracy"])
        m_integrator.Parameters().rtol = config["Initialize"]["Accuracy"].as<double>();

    return true;
}

void ProcessGroup::Optimize() {
    b_optimize = true;

    spdlog::info("Optimizing process group: Nucleus = {}, Nuclear Model = {}, Multiplicity = {}",
                 m_nucleus->ToString(), m_backend->GetNuclearModel()->GetName(),
                 m_processes[0].Info().Multiplicity());

    auto func = [&](const std::vector<FourVector> &mom, const double &wgt) {
        return SingleEvent(mom, wgt).Weight();
    };
    m_integrand.Function() = func;
    m_integrator.Optimize(m_integrand);

    spdlog::info("Calculating maximum weight");
    b_calc_weights = true;
    m_integrator.Parameters().ncalls = 100000;
    for(size_t i = 0; i < 3; ++i) m_integrator(m_integrand);
    b_calc_weights = false;
    b_optimize = false;

    // Store max weight and weight vector
    m_process_weights.resize(m_processes.size());
    for(size_t i = 0; i < m_processes.size(); ++i) {
        m_process_weights[i] = m_processes[i].MaxWeight();
        m_maxweight += m_process_weights[i];
        spdlog::info("Process xsec: {} ", m_processes[i].TotalCrossSection());
    }
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
    Event event = Event(m_nucleus.get(), mom, ps_wgt);

    spdlog::debug("Event Phase Space:");
    size_t idx = 0;
    for(const auto &momentum : event.Momentum()) {
        spdlog::debug("\t{}: {} (M2 = {})", ++idx, momentum, momentum.M2());
    }

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
