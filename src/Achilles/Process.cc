// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

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
#include "Plugins/Sherpa/SherpaInterface.hh"
#else
// Create dummy class
namespace achilles {
class SherpaInterface {};
} // namespace achilles
#endif

using achilles::Particle;
using achilles::Process;
using achilles::ProcessGroup;
using achilles::refParticles;
using achilles::FourVector;
using achilles::PID;
using achilles::ParticleStatus;

namespace achilles {
namespace {

using PIDs = std::vector<PID>;

// options[i] = the concrete PIDs allowed in slot i
std::vector<PIDs> CartesianProduct(const std::vector<PIDs> &options) {
    std::vector<PIDs> result{{}}; // one empty tuple
    for(const auto &opts : options) {
        std::vector<PIDs> next;
        next.reserve(result.size() * opts.size());
        for(const auto &partial : result)
            for(const auto &pid : opts) {
                auto extended = partial;
                extended.push_back(pid);
                next.push_back(std::move(extended));
            }
        result = std::move(next);
    }
    return result;
}

PIDs Options(PID pid) {
    ParticleInfo info(pid);
    if(!info.IsGroup()) return {pid};
    PIDs out;
    // IntID keeps the anti-particle sign; ID() would turn pi- into pi+
    for(const auto &m : info.Decompose()) out.push_back(PID{m.IntID()});
    return out;
}

std::vector<PIDs> ExpandSlots(const PIDs &slots) {
    std::vector<PIDs> options;
    for(auto pid : slots) options.push_back(Options(pid));
    return CartesianProduct(options);
}

std::vector<std::pair<PIDs, PIDs>> ExpandHadronic(const PIDs &in, const PIDs &out) {
    std::vector<PIDs> options;
    for(auto pid : in) options.push_back(Options(pid));
    for(auto pid : out) options.push_back(Options(pid));

    std::vector<std::pair<PIDs, PIDs>> expanded;
    for(const auto &combo : CartesianProduct(options)) {
        auto split = combo.begin() + static_cast<std::ptrdiff_t>(in.size());
        expanded.emplace_back(PIDs(combo.begin(), split), PIDs(split, combo.end()));
    }
    return expanded;
}

// Hadronic charge change in thirds: sum(final) - sum(initial)
int HadronicCharge3(const ProcessInfo &process) {
    int charge = 0;
    for(const auto &pid : process.m_hadronic.first) charge -= ParticleInfo(pid).IntCharge();
    for(const auto &pid : process.m_hadronic.second) charge += ParticleInfo(pid).IntCharge();
    return charge;
}

// Total charge change must vanish: ΔQ_lep + ΔQ_had = 0
bool ConservesCharge(const ProcessInfo &process) {
    return 3 * process.LeptonicCharge() + HadronicCharge3(process) == 0;
}

std::vector<ProcessInfo> ExpandProcessInfo(const ProcessInfo &process) {
    std::vector<ProcessInfo> result;
    const auto &[in, out] = process.m_hadronic;
    const auto spectators = ExpandSlots(process.m_spectator);

    for(auto &[new_in, new_out] : ExpandHadronic(in, out)) {
        ProcessInfo copy = process;
        copy.m_hadronic = {new_in, new_out};
        // TODO: Add MinimalFlips filter; MEC currently includes np->pn / pn->np for ΔQ = 0
        if(!ConservesCharge(copy)) continue;
        // Spectators do not change charge, so every combination is kept
        for(const auto &spectator : spectators) {
            copy.m_spectator = spectator;
            result.push_back(copy);
        }
    }
    return result;
}

} // namespace
} // namespace achilles

achilles::FourVector Process::ExtractQ(const Event &event) const {
    auto q = event.leptonIn().Momentum();
    for(const Particle& lepOut:event.leptonsOut())
		q-=lepOut.Momentum();
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

/*void ProcessGroup::SetupLeptons(Event &event, std::optional<size_t> process_idx) const {
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
}*/

void ProcessGroup::CrossSection(Event &event, std::optional<size_t> process_idx) {
    if(!process_idx) {
        double weight = 0;
        for(size_t i = 0; i < m_processes.size(); ++i) {
            double process_weight = m_backend->CrossSection(event, m_processes[i]);
            if(b_calc_weights) m_processes[i].AddWeight(process_weight);
            weight += process_weight;
        }
        event.Weight() = weight;
        if(b_calc_weights) m_xsec += weight;
    } else {
        Process& process = m_processes[process_idx.value()];
        double weight = m_backend->CrossSection(event, process);
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
        ProcessInfo process_info = process_node.as<ProcessInfo>();
        model->AllowedStates(process_info);
        const size_t multiplicity = process_info.Multiplicity();
        if(groups.find(multiplicity) == groups.end()) {
            groups.insert({multiplicity, ProcessGroup(beam, nucleus)});
        }
        ProcessGroup& group = groups.at(multiplicity);
        group.SetProcessInfo(process_info);
        group.BuildProcess(process_info, settings, model);
    }

    return groups;
}

void ProcessGroup::BuildProcess(const ProcessInfo &info, const Settings &settings,
                                NuclearModel *model) {
    const auto unweight_name = settings.GetAs<std::string>("Options/Unweighting/Name");

    const auto expanded_infos = ExpandProcessInfo(info);
    if(expanded_infos.empty())
        throw std::runtime_error(
            fmt::format("ProcessGroup: No charge-conserving states found for {}", info));

    for(const auto &expanded : expanded_infos) {
        try {
            auto unweighter =
                UnweighterFactory::Initialize(unweight_name, settings["Options/Unweighting"]);
            Process process(expanded, std::move(unweighter));
            process.SetID(model);
            AddProcess(std::move(process));
        } catch(std::out_of_range &e) {
            spdlog::error("Unweighter: Requested unweighter \"{}\", did you mean \"{}\"",
                          unweight_name,
                          achilles::GetSuggestion(UnweighterFactory::List(), unweight_name));
            exit(-1);
        }
    }
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

    // Build the phase space from the container-level process (container masses), so it does not
    // depend on which expanded process comes first. Groups built without one, or backends that
    // need concrete flavours, fall back to the first process.
    const ProcessInfo& channel_info =
        m_process_info && m_backend->ContainerChannels() ? *m_process_info : m_processes[0].Info();
    try {
        m_backend->SetupChannels(channel_info, m_beam, m_integrand, m_nucleus->ID());
    } catch(const InvalidChannel &) { return false; }

    UnitsEnum display_unit = UnitsEnum::nb;
    if(config.Exists("Main/DisplayUnit")) {
        display_unit = ToUnitEnum(config.GetAs<std::string>("Main/DisplayUnit"));
    }

    MultiChannelParams multichannel_params;
    if(config.Exists("Options/Initialize/Parameters"))
        multichannel_params = config.GetAs<MultiChannelParams>("Options/Initialize/Parameters");
    m_integrator = MultiChannel(m_integrand.NDims(), m_integrand.NChannels(), multichannel_params,
                                display_unit);
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

    auto func = [&](const Event& event, const double& wgt) {
		spdlog::trace("Cloning event...");
		Event clone=event; // because "event" is immutable. Is there a better way to do this?
		clone.Weight()=wgt;
		EventSetup(clone);
        return clone.Weight();
    };
    m_integrand.Function() = func;
    if(NeedsOptimization()) {
        spdlog::info(
            "Optimizing process group: Nucleus = {}, Nuclear Model = {}, Multiplicity = {}",
            m_nucleus->ToString(), m_backend->GetNuclearModel()->GetName(),
            m_processes[0].Info().Multiplicity());
		Event event(GroupProcessInfo(),m_nucleus);
        m_integrator.Optimize(m_integrand,event);
    }

    // Ensure that the integrator does not save these results into the summary
    m_integrator.UpdateSummary(false);
    if(m_maxweight == 0) {
        spdlog::info("Calculating maximum weight");
        b_calc_weights = true;
        m_integrator.Parameters().ncalls = 100000;

        for(size_t i = 0; i < 3; ++i) {
			// Create dummy point for training purposes
			Event event(GroupProcessInfo(),m_nucleus);
			m_integrator(m_integrand,event);
		}

        b_calc_weights = false;
        std::ofstream outFile("Results.txt");

        // Store max weight and weight vector
        m_process_weights.resize(m_processes.size());
        for(size_t i = 0; i < m_processes.size(); ++i) {
            m_process_weights[i] = m_processes[i].MaxWeight();
            m_maxweight += m_process_weights[i];
            spdlog::info("Process xsec: {} ", m_processes[i].TotalCrossSection());
            outFile << m_processes[i].TotalCrossSection() << "\n";
        }

        outFile.close();
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
	Event event(GroupProcessInfo(),m_nucleus);
	m_integrator.GeneratePoint(m_integrand,event);
	event.Weight()=m_integrator.GenerateWeight(m_integrand,event);
	EventSetup(event);
	return event;
}

void ProcessGroup::EventSetup(Event& event) {
    // Cut on leptons: NOTE: This assumes that all processes in the group have the same leptons
    if(!m_cuts.EvaluateCuts(event)) {
        // Ensure process weights are tracked correctly
        if(b_calc_weights) {
            for(Process& process : m_processes) process.AddWeight(0);
            m_xsec += 0;
        }
        event.Weight() = 0;
        return;
    }

    std::optional<size_t> process_opt = b_optimize ? std::nullopt : std::optional<size_t>(SelectProcess());
    CrossSection(event, process_opt);

    // If training the integrator or weight is zero, we can stop here
    if(b_optimize || event.Weight() == 0) return;

    // Otherwise, we need to fill the event with the selected process
    Process& process = m_processes[process_opt.value()];
    event.Flux() = m_beam->EvaluateFlux(process.Info().m_leptonic.first, event.leptonIn().Momentum());
    event.ProcessId() = process.ID();
	event.setProcessInfo(process.Info());
	// At this point, the event's PIDs have all been resolved;
	// we can now fill the rest of the event with anything PID-dependent,
	// such as the Remnant and Cascade
	SetupRemnant(event);
}

void ProcessGroup::SetupRemnant(Event& event) {
	EventVertex* targetVertex=event.uniqueVertex(EventVertex::StatusCode::target);
	Particle& init_nuc=event.particlesIn(targetVertex->Index())[0].get();
	if(!init_nuc.IsComposite())
		return;

	// Hadron and Remnant
	// hnr.first: Hadron to map to ev_hadron
	// hnr.second: Remnant after the split
	std::pair<Particle,Particle> hnr;
	Particle remnant=init_nuc;
	bool hasFoundPos=false;
	for(Particle& ev_hadron:event.particlesOut(targetVertex->Index())) {
		hnr=remnant.splitSingleRandom(ev_hadron.ID());
		remnant=hnr.second;

		// Splitting removes the hadron's momentum from the remnant by default,
		// but we're also changing that hadron's momentum (via mappers),
		// so we need to manually change the remnant's recoil.
		remnant.Momentum()+=hnr.first.Momentum()-ev_hadron.Momentum();
		ev_hadron.Position()=hnr.first.Position();

		if(!hasFoundPos) {
			event.Primary()->Position()=hnr.first.Position();
			targetVertex->Position()=hnr.first.Position();
			hasFoundPos=true;
		}
	}

	// Remaining bound particles are now tracked by the remnant;
	// composite nucleus no longer needed
	init_nuc.SubParts().clear();

	remnant.Status()=ParticleStatus::intermediate_residue;
	event.addParticle(remnant,targetVertex->Index(),Particle::OPEN_END);
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