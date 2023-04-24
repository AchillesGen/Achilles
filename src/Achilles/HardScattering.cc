#include <iostream>
#include <utility>

#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/HardScattering.hh"
#include "Achilles/HardScatteringFactory.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Random.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#include "METOOLS/Main/Spin_Structure.H"
#pragma GCC diagnostic pop
#endif

// Aliases for most common types
using achilles::HardScattering;
using achilles::Particles;

void HardScattering::SetProcess(const ProcessInfo &process) {
    spdlog::debug("Adding Process: {}", process);
    m_leptonicProcess = process;
#ifndef ACHILLES_SHERPA_INTERFACE
    m_current.Initialize(process);
    SMFormFactor = m_current.GetFormFactor();
#endif
}

std::vector<double> HardScattering::CrossSection(Event &event) const {
    // Calculate leptonic currents
    auto leptonCurrent = m_current.CalcCurrents(event.Momentum()[0], event.Momentum()[2]);

    // Calculate the hadronic currents
    // TODO: Clean this up and make generic for the nuclear model
    // TODO: Move this to initialization to remove check each time
    static std::vector<NuclearModel::FFInfoMap> ffInfo;
    if(ffInfo.empty()) {
        ffInfo.resize(3);
        for(const auto &current : leptonCurrent) {
#ifdef ACHILLES_SHERPA_INTERFACE
            ffInfo[0][current.first] = p_sherpa->FormFactors(PID::proton(), current.first);
            ffInfo[1][current.first] = p_sherpa->FormFactors(PID::neutron(), current.first);
            ffInfo[2][current.first] = p_sherpa->FormFactors(PID::carbon(), current.first);
#else
            // TODO: Define values somewhere
            ffInfo[0][current.first] = SMFormFactor.at({PID::proton(), current.first});
            ffInfo[1][current.first] = SMFormFactor.at({PID::neutron(), current.first});
            ffInfo[2][current.first] = SMFormFactor.at({PID::carbon(), current.first});
#endif
        }
    }
    auto hadronCurrent = m_nuclear->CalcCurrents({event.Momentum()[1]}, {event.Momentum()[3]},
                                                 event.Momentum()[0], ffInfo[0]);

    std::vector<double> amps2(hadronCurrent.size());
    const size_t nlep_spins = leptonCurrent.begin()->second.size();
    const size_t nhad_spins = m_nuclear->NSpins();

#ifdef ACHILLES_SHERPA_INTERFACE
    std::vector<METOOLS::Spin_Amplitudes> spin_amps;
    p_sherpa->FillAmplitudes(spin_amps);
    for(auto &amp : spin_amps)
        for(auto &elm : amp) elm = 0;
#else
    std::vector<std::vector<std::complex<double>>> spin_amps;
    spin_amps.resize(1);
    spin_amps[0].resize(nlep_spins * nhad_spins);
#endif

    for(size_t i = 0; i < nlep_spins; ++i) {
        // TODO: Fix this to be correct!!!
        size_t idx = ((i & ~1ul) << 2) + ((i & 1ul) << 1);
        idx = idx == 2 ? 10 : 2;
        // idx = 2
        for(size_t j = 0; j < nhad_spins; ++j) {
            double sign = 1.0;
            std::vector<std::complex<double>> amps(hadronCurrent.size());
            for(size_t mu = 0; mu < 4; ++mu) {
                for(const auto &lcurrent : leptonCurrent) {
                    auto boson = lcurrent.first;
                    if(hadronCurrent.find(boson) != hadronCurrent.end()) {
                        amps[0] += sign * lcurrent.second[i][mu] * hadronCurrent[boson][j][mu];
                        spin_amps[0][idx] +=
                            sign * lcurrent.second[i][mu] * hadronCurrent[boson][j][mu];
                    }
                }
                sign = -1.0;
            }
            for(size_t k = 0; k < hadronCurrent.size(); ++k) amps2[k] += std::norm(amps[k]);
        }
    }

#ifdef ACHILLES_SHERPA_INTERFACE
    for(const auto &amp : spin_amps) spdlog::trace("\n{}", amp);
    p_sherpa->FillAmplitudes(spin_amps);
#endif

    double spin_avg = 1;
    if(!ParticleInfo(m_leptonicProcess.m_leptonic.first).IsNeutrino()) spin_avg *= 2;
    if(m_nuclear->NSpins() > 1) spin_avg *= 2;

    // TODO: Correct this flux
    // double flux = 4*sqrt(pow(event.Momentum()[0]*event.Momentum()[1], 2)
    //                      -
    //                      event.Momentum()[0].M2()*event.Momentum()[1].M2());
    // TODO: Handle multiple hadron initial states
    double mass = ParticleInfo(m_leptonicProcess.m_hadronic.first[0]).Mass();
    double flux = 2 * event.Momentum()[1].E() * 2 * sqrt(event.Momentum()[0].P2() + mass * mass);
    static constexpr double to_nb = 1e6;
    std::vector<double> xsecs(hadronCurrent.size());
    for(size_t i = 0; i < hadronCurrent.size(); ++i) {
        xsecs[i] = amps2[i] * Constant::HBARC2 / spin_avg / flux * to_nb;
        spdlog::debug("Xsec[{}] = {}", i, xsecs[i]);
    }

    return xsecs;
}

bool HardScattering::FillEvent(Event &event, const std::vector<double> &xsecs) const {
    if(!m_nuclear->FillNucleus(event, xsecs)) return false;

    event.InitializeLeptons(m_leptonicProcess);
    event.InitializeHadrons(m_leptonicProcess);

    return true;
}
