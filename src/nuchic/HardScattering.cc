#include <iostream>
#include <utility>

#include "nuchic/HardScattering.hh"
#include "nuchic/Constants.hh"
#include "nuchic/NuclearModel.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/HardScatteringFactory.hh"
#include "nuchic/Event.hh"
#include "nuchic/Random.hh"

// Aliases for most common types
using nuchic::Particles;
using nuchic::HardScattering;

void HardScattering::SetProcess(const nuchic::Process_Info &process) {
    spdlog::debug("Adding Process: {}", process);
    m_leptonicProcess = process;
}

nuchic::Currents HardScattering::LeptonicCurrents(const std::vector<FourVector> &p,
                                                  const double &mu2) const {
    std::vector<std::array<double, 4>> mom(p.size());
    auto qVec = p[1];
    for(size_t i = 2; i < p.size()-1; ++i) {
        qVec -= p[i];
    }
    // auto rotMat = qVec.AlignZ();
    mom[0] = (-qVec/1_GeV).Momentum(); // .Rotate(rotMat).Momentum();
    mom[1] = (-p[1]/1_GeV).Momentum(); // .Rotate(rotMat).Momentum();
    for(size_t i = 2; i < p.size()-1; ++i) {
        mom[i] = (p[i]/1_GeV).Momentum(); // .Rotate(rotMat).Momentum();
    }
    mom.back() = {};

    // TODO: Clean up to handle the dummy hadronic states
    size_t idx = 0;
    std::vector<int> pids;
    pids.emplace_back(m_leptonicProcess.m_states.begin()->first[0]);
    spdlog::debug("PID: {}, Momentum: ({}, {}, {}, {})", pids[idx],
                  mom[idx][0], mom[idx][1], mom[idx][2], mom[idx][3]); 
    ++idx;
    for(const auto &pid : m_leptonicProcess.m_ids) {
        pids.emplace_back(pid);
        spdlog::debug("PID: {}, Momentum: ({}, {}, {}, {})", pids[idx],
                      mom[idx][0], mom[idx][1], mom[idx][2], mom[idx][3]); 
        ++idx;
    }
    pids.emplace_back(m_leptonicProcess.m_states.begin()->second[0]);
    spdlog::debug("PID: {}, Momentum: ({}, {}, {}, {})", pids[idx],
                  mom[idx][0], mom[idx][1], mom[idx][2], mom[idx][3]); 
    auto currents = p_sherpa -> Calc(pids, mom, mu2);

    for(auto &current : currents) { 
        spdlog::trace("Current for {}", current.first);
        for(size_t i = 0; i < current.second.size(); ++i) {
            for(size_t j = 0; j < current.second[0].size(); ++j) {
                if(mom.size() == 4) current.second[i][j] /= 1_GeV;
                else if(mom.size() == 6) current.second[i][j] /= pow(1_GeV, 3);
                spdlog::trace("Current[{}][{}] = {}", i, j, current.second[i][j]);
            }
        }
    }

    return currents;
}

std::vector<double> HardScattering::CrossSection(Event &event) const {
    // Calculate leptonic currents
    auto leptonCurrent = LeptonicCurrents(event.Momentum(), 100);

    // Calculate the hadronic currents
    // TODO: Clean this up and make generic for the nuclear model
    // TODO: Move this to initialization to remove check each time
    static std::vector<NuclearModel::FFInfoMap> ffInfo;
    if(ffInfo.empty()) {
        ffInfo.resize(3);
        for(const auto &current : leptonCurrent) {
            ffInfo[0][current.first] = p_sherpa -> FormFactors(PID::proton(), current.first);
            ffInfo[1][current.first] = p_sherpa -> FormFactors(PID::neutron(), current.first);
            ffInfo[2][current.first] = p_sherpa -> FormFactors(PID::carbon(), current.first);
        }
    }
    auto hadronCurrent = m_nuclear -> CalcCurrents(event, ffInfo);

    std::vector<double> amps2(hadronCurrent.size());
    const size_t nlep_spins = 1 << (event.Momentum().size() - 2);
    const size_t nhad_spins = m_nuclear -> NSpins();
    for(size_t i = 0; i  < nlep_spins; ++i) {
        for(size_t j = 0; j < nhad_spins; ++j) {
            double sign = 1.0;
            std::vector<std::complex<double>> amps(hadronCurrent.size());
            for(size_t mu = 0; mu < 4; ++mu) {
                for(const auto &lcurrent : leptonCurrent) {
                    auto boson = lcurrent.first;
                    for(size_t k = 0; k < hadronCurrent.size(); ++k) {
                        if(hadronCurrent[k].find(boson) != hadronCurrent[k].end())
                            amps[k] += sign*lcurrent.second[i][mu]*hadronCurrent[k][boson][j][mu];
                    }
                }
                sign = -1.0;
            }
            for(size_t k = 0; k < hadronCurrent.size(); ++k)
                amps2[k] += std::norm(amps[k]);
        }
    }

    double spin_avg = 1;
    if(!ParticleInfo(m_leptonicProcess.m_ids[0]).IsNeutrino()) spin_avg *= 2;
    if(m_nuclear -> NSpins() > 1) spin_avg *= 2;

    // TODO: Correct this flux
    // double flux = 4*sqrt(pow(event.Momentum()[0]*event.Momentum()[1], 2) 
    //                      - event.Momentum()[0].M2()*event.Momentum()[1].M2());
    double mass = ParticleInfo(m_leptonicProcess.m_states.begin()->first[0]).Mass();
    double flux = 2*event.Momentum()[1].E()*2*sqrt(event.Momentum()[0].P2() + mass*mass);
    static constexpr double to_nb = 1e6;
    std::vector<double> xsecs(hadronCurrent.size());
    for(size_t i = 0; i < hadronCurrent.size(); ++i) {
        xsecs[i] = amps2[i]*Constant::HBARC2/spin_avg/flux*to_nb;
        spdlog::debug("Xsec[{}] = {}", i, xsecs[i]);
    }

    return xsecs;
}

bool HardScattering::FillEvent(Event& event, const std::vector<double> &xsecs) const {
    if(!m_nuclear -> FillNucleus(event, xsecs))
        return false;
    
    event.InitializeLeptons(m_leptonicProcess);
    event.InitializeHadrons(m_leptonicProcess);

    return true;
}
