#include "Achilles/XSecBackend.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/NuclearModel.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/SherpaInterface.hh"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#include "METOOLS/Main/Spin_Structure.H"
#pragma GCC diagnostic pop
#endif

#include <stdexcept>

using achilles::XSecBackend;
using achilles::XSecBackendFactory;
using achilles::XSecBuilder;

void XSecBackend::ExtractMomentum(const Event &event, const ProcessInfo &process_info,
                                  FourVector &lep_in, std::vector<FourVector> &had_in,
                                  std::vector<FourVector> &lep_out,
                                  std::vector<FourVector> &had_out) const {
    static constexpr size_t lepton_in_end_idx = 1;
    const size_t hadron_end_idx = process_info.m_hadronic.first.size() + lepton_in_end_idx;
    const size_t lepton_end_idx = process_info.m_leptonic.second.size() + hadron_end_idx;
    auto momentum = event.Momentum();
    lep_in = momentum[0];
    had_in = std::vector<FourVector>(momentum.begin() + lepton_in_end_idx,
                                     momentum.begin() + static_cast<int>(hadron_end_idx));
    lep_out = std::vector<FourVector>(momentum.begin() + static_cast<int>(hadron_end_idx),
                                      momentum.begin() + static_cast<int>(lepton_end_idx));
    had_out = std::vector<FourVector>(momentum.begin() + static_cast<int>(lepton_end_idx),
                                      momentum.end());
}

double XSecBackend::FluxFactor(const FourVector &lep_in, const FourVector &had_in,
                               const ProcessInfo &process_info) const {
    double spin_avg = 1;
    if(!ParticleInfo(process_info.m_leptonic.first).IsNeutrino()) spin_avg *= 2;
    // TODO: Handle MEC correctly
    if(m_model->NSpins() > 1) spin_avg *= 2;

    // TODO: Correct this flux
    // double flux = 4*sqrt(pow(event.Momentum()[0]*event.Momentum()[1], 2)
    //                      -
    //                      event.Momentum()[0].M2()*event.Momentum()[1].M2());
    // TODO: Handle multiple hadron initial states
    double mass = ParticleInfo(process_info.m_hadronic.first[0]).Mass();
    double flux = 2 * lep_in.E() * 2 * sqrt(had_in.P2() + mass * mass);
    static constexpr double to_nb = 1e6;
    return Constant::HBARC2 / spin_avg / flux * to_nb;
}

achilles::DefaultBackend::DefaultBackend() {}

double achilles::DefaultBackend::CrossSection(const Event &event, const Process &process) const {
    const auto &process_info = process.Info();
    FourVector lepton_momentum_in;
    std::vector<FourVector> hadron_momentum_in, hadron_momentum_out, lepton_momentum_out;
    ExtractMomentum(event, process_info, lepton_momentum_in, hadron_momentum_in,
                    lepton_momentum_out, hadron_momentum_out);

    std::pair<PID, PID> leptons{process_info.m_leptonic.first, process_info.m_leptonic.second[0]};
    auto lcurrent_calculator = m_currents.at(leptons);

    auto lepton_current =
        lcurrent_calculator.CalcCurrents(lepton_momentum_in, lepton_momentum_out[0]);

    // TODO: Handle the case for MEC, RES, and DIS
    NuclearModel::FFInfoMap ff_info;
    for(const auto &boson : lepton_current) {
        ff_info[boson.first] = form_factors.at({process_info.m_hadronic.first[0], boson.first});
    }
    auto q = lepton_momentum_in - lepton_momentum_out[0];
    auto hadron_current =
        m_model->CalcCurrents(hadron_momentum_in, hadron_momentum_out, q, ff_info);

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

    return amps2 * FluxFactor(lepton_momentum_in, hadron_momentum_in[0], process_info);
}

void achilles::DefaultBackend::AddProcess(Process &process) {
    const auto &process_info = process.Info();
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
    if(form_factors.size() == 0) form_factors = current.GetFormFactor();
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

double achilles::BSMBackend::CrossSection(const Event &event, const Process &process) const {
    const auto &process_info = process.Info();
    FourVector lepton_momentum_in;
    std::vector<FourVector> hadron_momentum_in, hadron_momentum_out, lepton_momentum_out;
    ExtractMomentum(event, process_info, lepton_momentum_in, hadron_momentum_in,
                    lepton_momentum_out, hadron_momentum_out);
    auto lepton_current = CalcLeptonCurrents(event.Momentum(), process_info);

    // TODO: Handle the case for MEC, RES, and DIS
    NuclearModel::FFInfoMap ff_info;
    for(const auto &boson : lepton_current) {
        ff_info[boson.first] = p_sherpa->FormFactors(process_info.m_hadronic.first[0], boson.first);
    }
    auto q = lepton_momentum_in;
    for(const auto &mom : lepton_momentum_out) q -= mom;
    auto hadron_current =
        m_model->CalcCurrents(hadron_momentum_in, hadron_momentum_out, q, ff_info);

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

    return amps2 * FluxFactor(lepton_momentum_in, hadron_momentum_in[0], process_info);
}

achilles::Currents achilles::BSMBackend::CalcLeptonCurrents(const std::vector<FourVector> &p,
                                                            const ProcessInfo &info) const {
    // TODO: Move adapter code into Sherpa interface code
    std::vector<std::array<double, 4>> mom(p.size());
    std::vector<int> pids;
    for(const auto &elm : info.m_mom_map) {
        pids.push_back(static_cast<int>(elm.second));
        mom[elm.first] = (p[elm.first] / 1_GeV).Momentum();
        spdlog::debug("PID: {}, Momentum: ({}, {}, {}, {})", pids.back(), mom[elm.first][0],
                      mom[elm.first][1], mom[elm.first][2], mom[elm.first][3]);
    }
    // TODO: Figure out if we want to have a scale dependence (Maybe for DIS??)
    static constexpr double mu2 = 100;
    auto currents = p_sherpa->Calc(pids, mom, mu2);

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

achilles::SherpaBackend::SherpaBackend() {}

void achilles::SherpaBackend::AddProcess(Process &process) {
    auto &process_info = process.Info();
    spdlog::debug("Initializing leptonic currents");
    if(!p_sherpa->InitializeProcess(process_info)) {
        spdlog::error("Cannot initialize hard process");
        exit(1);
    }
    process_info.m_mom_map = p_sherpa->MomentumMap(process_info.Ids());
}

double achilles::SherpaBackend::CrossSection(const Event &, const Process &) const {
    throw std::runtime_error("SherpaBackend: Currently not implemented");
    return 0;
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

XSecBuilder &XSecBuilder::AddNuclearModel(std::shared_ptr<NuclearModel> model) {
    m_backend->SetNuclearModel(model);
    return *this;
}

std::unique_ptr<XSecBackend> XSecBuilder::build() {
    if(!m_backend->Validate()) throw std::runtime_error("Backend is not properly configured!");
    return std::move(m_backend);
}
