#include "Achilles/XSecBackend.hh"
#include "Achilles/Channels.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Process.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/SherpaInterface.hh"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#include "METOOLS/Main/Spin_Structure.H"
#include "plugins/Sherpa/FormFactors.hh"
#pragma GCC diagnostic pop
#endif

#include <stdexcept>

using achilles::XSecBackend;
using achilles::XSecBackendFactory;
using achilles::XSecBuilder;

double XSecBackend::SpinAvg(const ProcessInfo &process_info) const {
    double spin_avg = 1;
    if(!ParticleInfo(process_info.m_leptonic.first).IsNeutrino()) spin_avg *= 2;
    // TODO: Handle MEC correctly
    if(m_model->NSpins() > 1) spin_avg *= 2;

    return 1.0 / spin_avg;
}

double XSecBackend::FluxFactor(const FourVector &lep_in, const FourVector &had_in,
                               const ProcessInfo &process_info) const {
    // TODO: Correct this flux
    // double flux = 4*sqrt(pow(event.Momentum()[0]*event.Momentum()[1], 2)
    //                      -
    //                      event.Momentum()[0].M2()*event.Momentum()[1].M2());
    // TODO: Handle multiple hadron initial states
    double mass = ParticleInfo(process_info.m_hadronic.first[0]).Mass();
    double flux = 2 * lep_in.E() * 2 * sqrt(had_in.P2() + mass * mass);
    static constexpr double to_nb = 1e6;
    return Constant::HBARC2 / flux * to_nb;
}

double XSecBackend::InitialStateFactor(size_t nprotons, size_t nneutrons,
                                       const std::vector<Particle> &p_in, const std::vector<Particle> &p_spect) const {
    auto initial_wgt = m_model->InitialStateWeight(p_in, p_spect, nprotons, nneutrons);
    return initial_wgt;
}

achilles::DefaultBackend::DefaultBackend() {}

double achilles::DefaultBackend::CrossSection(const Event &event_in, const Process &process) const {
    auto event = event_in;

    for(const auto &part : event.Momentum()) {
        spdlog::debug("Momentum: ({}, {}, {}, {})", part[0], part[1], part[2], part[3]);
    }
    m_model->TransformFrame(event, process, true);

    for(const auto &part : event.Momentum()) {
        spdlog::debug("Momentum: ({}, {}, {}, {})", part[0], part[1], part[2], part[3]);
    }

    const auto &process_info = process.Info();
    Particle lepton_in;
    std::vector<Particle> hadron_in, hadron_out, lepton_out, spect;
    process.ExtractParticles(event, lepton_in, hadron_in, lepton_out, hadron_out, spect);

    std::pair<PID, PID> leptons{process_info.m_leptonic.first, process_info.m_leptonic.second[0]};
    auto lcurrent_calculator = m_currents.at(leptons);

    auto lepton_current =
        lcurrent_calculator.CalcCurrents(lepton_in.Momentum(), lepton_out[0].Momentum());

    // TODO: Handle the case for MEC, RES, and DIS
    NuclearModel::FFInfoMap ff_info;
    for(const auto &boson : lepton_current) {
        ff_info[boson.first] = form_factors.at({process_info.m_hadronic.first[0], boson.first});
    }
    auto q = lepton_in.Momentum() - lepton_out[0].Momentum();
    auto hadron_current = m_model->CalcCurrents(hadron_in, hadron_out, spect, q, ff_info);

    double amps2 = 0;

    for(const auto &lcurrent : lepton_current) {
        auto boson = lcurrent.first;
        if(hadron_current.find(boson) == hadron_current.end()) continue;
        const auto &hcurrent = hadron_current[boson];
        for(const auto &lcurrent_spin : lcurrent.second) {
            for(const auto &hcurrent_spin : hcurrent) {
                amps2 += std::norm(lcurrent_spin * hcurrent_spin);
            }
        }
    }
    if(std::isnan(amps2)) amps2 = 0;

    auto flux = FluxFactor(lepton_in.Momentum(), hadron_in[0].Momentum(), process_info);
    size_t nprotons = event.CurrentNucleus()->NProtons();
    size_t nneutrons = event.CurrentNucleus()->NNeutrons();
    auto initial_wgt = InitialStateFactor(nprotons, nneutrons, hadron_in, spect);
    spdlog::debug("flux = {}, initial_wgt = {}, amps2 = {}", flux, initial_wgt,
                  amps2 * SpinAvg(process_info));
    double xsec = amps2 * flux * initial_wgt * SpinAvg(process_info) * event.Weight();
    if(std::isnan(xsec)) {
        spdlog::warn("Got nan for xsec, setting to 0");
        xsec = 0;
    }
    return xsec;
}

void achilles::DefaultBackend::AddProcess(Process &process) {
    const auto &process_info = process.Info();
    spdlog::debug("Adding process: {}", process.Info());
    if(process_info.m_leptonic.second.size() != 1)
        throw std::runtime_error("Achilles::DefaultBackend: Can only handle "
                                 "single lepton final states");

    std::pair<PID, PID> leptons{process_info.m_leptonic.first, process_info.m_leptonic.second[0]};
    if(m_currents.find(leptons) != m_currents.end()) {
        spdlog::debug("Achilles::DefaultBackend: Leptonic current for {} -> {} "
                      "already created",
                      leptons.first, leptons.second);
        return;
    }

    LeptonicCurrent current;
    current.Initialize(process_info);

    // TODO: Clean up how the form factors are loaded. Make part of the backend?
    if(form_factors.size() == 0)
        form_factors = current.GetFormFactor();
    else {
        auto tmp = current.GetFormFactor();
        form_factors.insert(tmp.begin(), tmp.end());
    }
    m_currents[leptons] = current;
}

void achilles::DefaultBackend::SetupChannels(const ProcessInfo &process_info,
                                             std::shared_ptr<Beam> beam,
                                             Integrand<FourVector> &integrand, PID nuc_id) {
    auto masses = process_info.Masses();
    auto multiplicity = process_info.FinalStateMultiplicity();
    if(multiplicity > 3) {
        const std::string error =
            fmt::format("Leptonic Tensor can only handle n->2 and n->3 processes without "
                        "BSM being enabled. "
                        "Got a n->{} process",
                        multiplicity);
        throw std::runtime_error(error);
    }

    if(multiplicity == 2) {
        Channel<FourVector> channel0 =
            BuildChannel<TwoBodyMapper>(m_model.get(), 2, 2, process_info.m_spectator.size(), beam, masses, nuc_id);
        integrand.AddChannel(std::move(channel0));
    } else {
        Channel<FourVector> channel0 =
            BuildChannel<ThreeBodyMapper>(m_model.get(), 3, 2, process_info.m_spectator.size(), beam, masses, nuc_id);
        integrand.AddChannel(std::move(channel0));
    }
}

#ifdef ACHILLES_SHERPA_INTERFACE
achilles::BSMBackend::BSMBackend() {}

void achilles::BSMBackend::AddProcess(Process &process) {
    auto &process_info = process.Info();
    spdlog::debug("Initializing leptonic currents");
    if(!p_sherpa->InitializeProcess(process_info)) {
        spdlog::error("Cannot initialize hard process");
        exit(1);
    }
    process_info.m_mom_map = p_sherpa->MomentumMap(process_info.Ids());
}

double achilles::BSMBackend::CrossSection(const Event &event_in, const Process &process) const {
    auto event = event_in;
    m_model->TransformFrame(event, process, true);

    const auto &process_info = process.Info();
    Particle lepton_in;
    std::vector<Particle> hadron_in, hadron_out, lepton_out, spect;
    process.ExtractParticles(event, lepton_in, hadron_in, lepton_out, hadron_out, spect);
    auto lepton_current = CalcLeptonCurrents(event.Momentum(), process_info);

    // TODO: Handle the case for MEC, RES, and DIS
    NuclearModel::FFInfoMap ff_info;
    for(const auto &boson : lepton_current) {
        ff_info[boson.first] = p_sherpa->FormFactors(process_info.m_hadronic.first[0], boson.first);
    }
    auto q = lepton_in.Momentum();
    for(const auto &part : lepton_out) q -= part.Momentum();
    auto hadron_current = m_model->CalcCurrents(hadron_in, hadron_out, spect, q, ff_info);

    // Setup handling of spin decays
    const size_t nlep_spins = lepton_current.begin()->second.size();
    const size_t nhad_spins = m_model->NSpins();
    std::vector<METOOLS::Spin_Amplitudes> spin_amps;
    p_sherpa->FillAmplitudes(spin_amps);
    for(auto &amp : spin_amps)
        for(auto &elm : amp) elm = 0;

    double amps2 = 0;
    for(const auto &lcurrent : lepton_current) {
        auto boson = lcurrent.first;
        if(hadron_current.find(boson) == hadron_current.end()) continue;
        const auto &hcurrent = hadron_current[boson];
        for(size_t i = 0; i < nlep_spins; ++i) {
            // FIXME: this to be correct!!!
            size_t idx = ((i & ~1ul) << 2) + ((i & 1ul) << 1);
            idx = idx == 2 ? 10 : 2;
            for(size_t j = 0; j < nhad_spins; ++j) {
                amps2 += std::norm(lcurrent.second[i] * hcurrent[j]);
                spin_amps[0][idx] += lcurrent.second[i] * hcurrent[j];
            }
        }
    }
    p_sherpa->FillAmplitudes(spin_amps);

    auto flux = FluxFactor(lepton_in.Momentum(), hadron_in[0].Momentum(), process_info);
    size_t nprotons = event.CurrentNucleus()->NProtons();
    size_t nneutrons = event.CurrentNucleus()->NNeutrons();
    auto initial_wgt = InitialStateFactor(nprotons, nneutrons, hadron_in, spect);
    return amps2 * flux * initial_wgt * SpinAvg(process_info) * event.Weight();
}

achilles::Currents achilles::BSMBackend::CalcLeptonCurrents(const std::vector<FourVector> &p,
                                                            const ProcessInfo &info) const {
    // TODO: Move adapter code into Sherpa interface code
    std::vector<std::array<double, 4>> mom(p.size());
    std::vector<int> pids;
    spdlog::debug("mom map: {}", info.m_mom_map.size());
    for(const auto &elm : info.m_mom_map) {
        pids.push_back(static_cast<int>(elm.second));
        mom[elm.first] = (p[elm.first] / 1_GeV).Momentum();
        spdlog::debug("PID: {}, Momentum: ({}, {}, {}, {})", pids.back(), mom[elm.first][0],
                      mom[elm.first][1], mom[elm.first][2], mom[elm.first][3]);
    }
    // TODO: Figure out if we want to have a scale dependence (Maybe for DIS??)
    static constexpr double mu2 = 100;
    auto currents = p_sherpa->CalcCurrent(pids, mom, mu2);

    for(auto &current : currents) {
        spdlog::trace("Current for {}", current.first);
        for(size_t i = 0; i < current.second.size(); ++i) {
            for(size_t j = 0; j < current.second[0].size(); ++j) {
                current.second[i][j] /= pow(1_GeV, static_cast<double>(mom.size()) - 3);
                spdlog::trace("Current[{}][{}] = {}", i, j, current.second[i][j]);
            }
        }
    }
    return currents;
}

void achilles::BSMBackend::SetupChannels(const ProcessInfo &process_info,
                                         std::shared_ptr<Beam> beam,
                                         Integrand<FourVector> &integrand, PID nuc_id) {
    auto masses = process_info.Masses();
    auto channels = p_sherpa->GenerateChannels(process_info.Ids());
    size_t count = 0;
    for(auto &chan : channels) {
        Channel<FourVector> channel =
            achilles::BuildGenChannel(m_model.get(), process_info.m_leptonic.second.size() + 1, 2,
                                      beam, std::move(chan), masses, nuc_id);
        integrand.AddChannel(std::move(channel));
        spdlog::info("Adding Channel{}", count++);
    }
}

void achilles::BSMBackend::SetOptions(const YAML::Node &options) {
    auto config = YAML::LoadFile(options["FormFactorFile"].as<std::string>());
    const auto vectorFF = config["vector"].as<std::string>();
    const auto axialFF = config["axial"].as<std::string>();
    const auto coherentFF = config["coherent"].as<std::string>();
    auto ff = FormFactorBuilder()
                  .Vector(vectorFF, config[vectorFF])
                  .AxialVector(axialFF, config[axialFF])
                  .Coherent(coherentFF, config[coherentFF])
                  .build();
    FormFactorInterface::SetFormFactor(std::move(ff));
}

achilles::SherpaBackend::SherpaBackend() {}

void achilles::SherpaBackend::AddProcess(Process &process) {
    auto &process_info = process.Info();
    spdlog::debug("Initializing leptonic currents");
    if(!p_sherpa->InitializeProcess(process_info)) {
        spdlog::error("Cannot initialize hard process");
        exit(1);
    }
}

double achilles::SherpaBackend::CrossSection(const Event &event, const Process &process) const {
    auto p = event.Momentum();
    auto info = process.Info();
    // TODO: Move adapter code into Sherpa interface code
    std::vector<std::array<double, 4>> mom(p.size());
    std::vector<long> pids = process.Info().Ids();
    for(size_t i = 0; i < event.Momentum().size(); ++i) { mom[i] = (p[i] / 1_GeV).Momentum(); }
    // TODO: Figure out if we want to have a scale dependence (Maybe for DIS??)
    static constexpr double mu2 = 100;
    auto amps2 = p_sherpa->CalcDifferential(pids, mom, mu2);
    spdlog::trace("|M|^2 = {}", amps2);
    if(std::isnan(amps2)) amps2 = 0;

    Particle lepton_in;
    std::vector<Particle> hadron_in, hadron_out, lepton_out, spect;
    process.ExtractParticles(event, lepton_in, hadron_in, lepton_out, hadron_out, spect);
    auto flux = FluxFactor(lepton_in.Momentum(), hadron_in[0].Momentum(), info);
    size_t nprotons = event.CurrentNucleus()->NProtons();
    size_t nneutrons = event.CurrentNucleus()->NNeutrons();
    auto initial_wgt = InitialStateFactor(nprotons, nneutrons, hadron_in, spect);
    spdlog::debug("flux = {}, initial_wgt = {}, amps2 = {}", flux, initial_wgt, amps2);
    return amps2 * flux * initial_wgt * event.Weight();
}

void achilles::SherpaBackend::SetupChannels(const ProcessInfo &process_info,
                                            std::shared_ptr<Beam> beam,
                                            Integrand<FourVector> &integrand, PID nuc_id) {
    auto masses = process_info.Masses();
    auto channels = p_sherpa->GenerateChannels(process_info.Ids());
    size_t count = 0;
    for(auto &chan : channels) {
        Channel<FourVector> channel =
            achilles::BuildGenChannel(m_model.get(), process_info.m_leptonic.second.size() + 1, 2,
                                      beam, std::move(chan), masses, nuc_id);
        integrand.AddChannel(std::move(channel));
        spdlog::info("Adding Channel{}", count++);
    }
}

void achilles::SherpaBackend::SetOptions(const YAML::Node &options) {
    auto config = YAML::LoadFile(options["FormFactorFile"].as<std::string>());
    const auto vectorFF = config["vector"].as<std::string>();
    const auto axialFF = config["axial"].as<std::string>();
    const auto coherentFF = config["coherent"].as<std::string>();
    auto ff = FormFactorBuilder()
                  .Vector(vectorFF, config[vectorFF])
                  .AxialVector(axialFF, config[axialFF])
                  .Coherent(coherentFF, config[coherentFF])
                  .build();
    FormFactorInterface::SetFormFactor(std::move(ff));
}
#endif

XSecBuilder::XSecBuilder(const std::string &name)
    : m_backend{XSecBackendFactory::Initialize(name)} {}

XSecBuilder &XSecBuilder::AddOptions(const YAML::Node &node) {
    m_backend->SetOptions(node);
    return *this;
}

XSecBuilder &XSecBuilder::AddSherpa(achilles::SherpaInterface *sherpa) {
    m_backend->SetSherpa(sherpa);
    return *this;
}

XSecBuilder &XSecBuilder::AddProcess(Process &process) {
    m_backend->AddProcess(process);
    return *this;
}

XSecBuilder &XSecBuilder::AddNuclearModel(std::unique_ptr<NuclearModel> model) {
    m_backend->AddNuclearModel(std::move(model));
    return *this;
}

std::unique_ptr<XSecBackend> XSecBuilder::build() {
    if(!m_backend->Validate()) throw std::runtime_error("Backend is not properly configured!");
    return std::move(m_backend);
}
