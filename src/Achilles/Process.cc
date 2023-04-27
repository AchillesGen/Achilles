#include "Achilles/Process.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Channels.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/SherpaInterface.hh"
#else
#include "Achilles/DummySherpa.hh"
#endif

using achilles::Process;
using achilles::ProcessGroup;

void Process::SelectInitialState(Event &event) const {
    FourVector lep_in;
    std::vector<FourVector> had_in, lep_out, had_out;
    ExtractMomentum(event, lep_in, had_in, lep_out, had_out);
    std::vector<Particle> leptons, hadrons;

    // Setup leptons
    leptons.emplace_back(m_info.m_leptonic.first, lep_in);
    leptons.back().Status() = ParticleStatus::initial_state;
    if(lep_out.size() != m_info.m_leptonic.second.size())
        throw std::runtime_error(
            "Process: Number of lepton momenta does not match number of leptons in event");
    for(size_t i = 0; i < lep_out.size(); ++i) {
        leptons.emplace_back(m_info.m_leptonic.second[i], lep_out[i]);
        leptons.back().Status() = ParticleStatus::final_state;
    }

    // Setup hadrons
    if(had_in.size() != m_info.m_hadronic.first.size() ||
       had_out.size() != m_info.m_hadronic.second.size())
        throw std::runtime_error(
            "Process: Number of hadrons momenta does not match number of hadrons in event");

    // Handle coherent scattering as a special case
    if(ParticleInfo(m_info.m_hadronic.first[0]).IsNucleus()) {
        event.CurrentNucleus()->Nucleons().clear();
        const auto id_in = m_info.m_hadronic.first[0];
        const auto id_out = m_info.m_hadronic.second[0];

        Particle initial(id_in, had_in[0]);
        initial.Status() = ParticleStatus::initial_state;
        event.CurrentNucleus()->Nucleons().push_back(initial);

        Particle final(id_out, had_out[0]);
        final.Status() = ParticleStatus::final_state;
        event.CurrentNucleus()->Nucleons().push_back(final);
        return;
    }

    // TODO: Handle position dependent probabilities
    auto protons = event.CurrentNucleus()->ProtonsIDs();
    auto neutrons = event.CurrentNucleus()->NeutronsIDs();
    Random::Instance().Shuffle(protons);
    Random::Instance().Shuffle(neutrons);
    std::vector<size_t> indices(had_in.size());
    for(size_t i = 0; i < had_in.size(); ++i) {
        if(m_info.m_hadronic.first[i] == PID::proton()) {
            indices[i] = protons.back();
            protons.pop_back();
        } else {
            indices[i] = neutrons.back();
            neutrons.pop_back();
        }

        auto &initial = event.CurrentNucleus()->Nucleons()[indices[i]];
        initial.Momentum() = had_in[i];
        initial.Status() = ParticleStatus::initial_state;
    }

    // TODO: Propagate deltas away from interaction point
    // TODO: Handle something like MEC+pion???
    ThreeVector position;
    size_t cur_idx = 0;
    for(size_t i = 0; i < had_out.size(); ++i) {
        Particle part;
        if(ParticleInfo(m_info.m_hadronic.second[i]).IsBaryon()) {
            part = event.CurrentNucleus()->Nucleons()[indices[cur_idx++]];
        } else {
            part = Particle(m_info.m_hadronic.second[i]);
        }
        part.Status() = ParticleStatus::propagating;
        part.Momentum() = had_out[i];
        position = part.Position();
        event.CurrentNucleus()->Nucleons().push_back(part);
    }
}

void Process::ExtractMomentum(const Event &event, FourVector &lep_in,
                              std::vector<FourVector> &had_in, std::vector<FourVector> &lep_out,
                              std::vector<FourVector> &had_out) const {
    static constexpr size_t lepton_in_end_idx = 1;
    const size_t hadron_end_idx = m_info.m_hadronic.first.size() + lepton_in_end_idx;
    const size_t lepton_end_idx = m_info.m_leptonic.second.size() + hadron_end_idx;
    auto momentum = event.Momentum();
    lep_in = momentum[0];
    had_in = std::vector<FourVector>(momentum.begin() + lepton_in_end_idx,
                                     momentum.begin() + static_cast<int>(hadron_end_idx));
    lep_out = std::vector<FourVector>(momentum.begin() + static_cast<int>(hadron_end_idx),
                                      momentum.begin() + static_cast<int>(lepton_end_idx));
    had_out = std::vector<FourVector>(momentum.begin() + static_cast<int>(lepton_end_idx),
                                      momentum.end());
}

double ProcessGroup::CrossSection(Event &event, std::optional<size_t> process_idx) {
    if(!process_idx) {
        double weight = 0;
        for(size_t i = 0; i < m_processes.size(); ++i) {
            auto process_weight = m_backend->CrossSection(event, m_processes[i]);
            m_processes[i].AddWeight(process_weight);
            weight += process_weight;
        }
        return weight;
    } else {
        auto &process = m_processes[process_idx.value()];
        auto weight = m_backend->CrossSection(event, process);
        weight = process.Unweight(weight);
        return weight;
    }
}

size_t ProcessGroup::SelectProcess() const {
    return Random::Instance().SelectIndex(m_process_weights);
}

std::map<size_t, ProcessGroup>
ProcessGroup::ConstructProcessGroups(const YAML::Node &node, std::shared_ptr<NuclearModel> model,
                                     std::shared_ptr<Beam> beam, std::shared_ptr<Nucleus> nucleus) {
    std::map<size_t, ProcessGroup> groups;
    for(const auto &process_node : node["Processes"]) {
        auto process_info = process_node.as<ProcessInfo>();
        model->AllowedStates(process_info);
        const auto unweight_name = node["Unweighting"]["Name"].as<std::string>();
        auto unweighter = UnweighterFactory::Initialize(unweight_name, node["Unweighting"]);
        Process process(process_info, std::move(unweighter));
        const auto multiplicity = process_info.Multiplicity();
        if(groups.find(multiplicity) == groups.end()) {
            groups[multiplicity] = ProcessGroup(model, beam, nucleus);
        }
        groups[multiplicity].AddProcess(std::move(process));
    }

    return groups;
}

void ProcessGroup::SetupBackend(const YAML::Node &node, SherpaInterface *sherpa) {
    auto backend_name = node["Backend"]["Name"].as<std::string>();
    auto &backend_builder = XSecBuilder(backend_name)
                                .AddOptions(node["Backend"]["Options"])
                                .AddNuclearModel(m_nuc_model)
                                .AddSherpa(sherpa);
    for(const auto &process : m_processes) {
        backend_builder = std::move(backend_builder.AddProcess(process));
    }
    m_backend = backend_builder.build();
}

void ProcessGroup::SetupIntegration(const YAML::Node &config) {
    if(m_processes.size() == 0)
        throw std::runtime_error("ProcessGroup: Number of processes found is zero!");

    m_backend->SetupChannels(m_processes[0].Info(), m_beam, m_integrand)

        // TODO: Fix scaling to be consistent with Chili paper
        m_integrator = MultiChannel(m_integrand.NDims(), m_integrand.NChannels(), {1000, 2});
    if(config["Initialize"]["Accuracy"])
        m_integrator.Parameters().rtol = config["Initialize"]["Accuracy"].as<double>();
}

void ProcessGroup::Optimize() {
    b_optimize = true;

    auto func = [&](const std::vector<FourVector> &mom, const double &wgt) {
        return SingleEvent(mom, wgt).Weight();
    };
    m_integrand.Function() = func;
    m_integrator.Optimize(m_integrand);
    m_integrator.Summary();

    // Store max weight and weight vector
    for(size_t i = 0; i < m_processes.size(); ++i) {
        m_process_weights[i] = m_processes[i].MaxWeight();
        m_maxweight += m_process_weights[i];
    }
    for(size_t i = 0; i < m_processes.size(); ++i) { m_process_weights[i] /= m_maxweight; }
    b_optimize = false;
}

void ProcessGroup::Summary() const {}

achilles::Event ProcessGroup::GenerateEvent() {
    const auto mom = m_integrator.GeneratePoint(m_integrand);
    const auto ps_wgt = m_integrator.GenerateWeight(m_integrand, mom);
    return SingleEvent(mom, ps_wgt);
}

achilles::Event ProcessGroup::SingleEvent(const std::vector<FourVector> &mom, double ps_wgt) {
    Event event(m_nucleus, mom, ps_wgt);

    spdlog::debug("Event Phase Space:");
    size_t idx = 0;
    for(const auto &momentum : event.Momentum()) {
        spdlog::debug("\t{}: {} (M2 = {})", ++idx, momentum, momentum.M2());
    }

    auto process_opt = b_optimize ? std::nullopt : std::optional<size_t>(SelectProcess());
    event.Weight() *= CrossSection(event, process_opt);

    // If training the integrator, we can stop here
    if(b_optimize) return event;

    // Otherwise, we need to fill the event with the selected process
    auto &process = m_processes[process_opt.value()];
    event.Flux() = m_beam->EvaluateFlux(process.Info().m_leptonic.first, event.Momentum()[0]);
    process.SelectInitialState(event);

    spdlog::trace("Leptons:");
    idx = 0;
    for(const auto &particle : event.Leptons()) { spdlog::trace("\t{}: {}", ++idx, particle); }

    spdlog::trace("Hadrons:");
    idx = 0;
    for(const auto &particle : event.Hadrons()) { spdlog::trace("\t{}: {}", ++idx, particle); }

    event.CalcWeight();
    spdlog::trace("Weight: {}", event.Weight());

    return event;
}
