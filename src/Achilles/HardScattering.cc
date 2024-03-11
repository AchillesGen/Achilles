#include <iostream>
#include <utility>

#include "Achilles/HardScattering.hh"
#include "Achilles/Constants.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/HardScatteringFactory.hh"
#include "Achilles/Event.hh"
#include "Achilles/Random.hh"

// Aliases for most common types
using achilles::Particles;
using achilles::HardScattering;
using achilles::LeptonicCurrent;

void LeptonicCurrent::Initialize(const Process_Info &process) {
    using namespace achilles::Constant;
    const std::complex<double> i(0, 1);
    // Determine process
    bool init_neutrino = ParticleInfo(process.ids[0]).IsNeutrino();
    bool neutral_current = NeutralCurrent(process.ids[0], process.ids[1]);
    bool charged_current = ChargedCurrent(init_neutrino, process.ids[0], process.ids[1]);
    if(!neutral_current && !charged_current) 
        throw std::runtime_error("HardScattering: Invalid process");

    // TODO: Define couplings correctly
    if(charged_current) {
        pid = init_neutrino ? (process.ids[0].AsInt() < 0 ? 24 : -24)
                            : (process.ids[0].AsInt() < 0 ? -24 : 24);
        coupl_right = 0;
        coupl_left = ee*i/(sw*sqrt(2));
        mass = Constant::MW;
        width = Constant::GAMW;
    } else if(neutral_current) {
        if(init_neutrino) {
            coupl_left = (cw*ee*i)/(2*sw)+(ee*i*sw)/(2*cw);
            coupl_right = 0;
            pid = 23;
            mass = Constant::MZ;
            width = Constant::GAMZ;
        } else {
            coupl_right = -ee*i;
            coupl_left = coupl_right;
            pid = 22;
        }
    }
    anti = process.ids[0].AsInt() < 0;
}

bool LeptonicCurrent::NeutralCurrent(achilles::PID initial, achilles::PID final) const {
    return initial == final;
}

bool LeptonicCurrent::ChargedCurrent(bool neutrino, achilles::PID initial, achilles::PID final) const {
    return std::abs(initial.AsInt()) - (2*neutrino - 1) == std::abs(final.AsInt());
}

achilles::FFDictionary LeptonicCurrent::GetFormFactor() {
    FFDictionary results;
    static constexpr std::complex<double> i(0, 1);
    using namespace achilles::Constant;
    // TODO: Double check form factors
    if(pid == 24) {
        const std::complex<double> coupl = ee*i/(sw*sqrt(2)*2);
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                         {FormFactorInfo::Type::F1n, -coupl},
                                         {FormFactorInfo::Type::F2p, coupl},
                                         {FormFactorInfo::Type::F2n, -coupl},
                                         {FormFactorInfo::Type::FA, coupl}};
        results[{PID::neutron(), pid}] = {};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == -24) {
        const std::complex<double> coupl = ee*i/(sw*sqrt(2)*2);
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                          {FormFactorInfo::Type::F1n, -coupl},
                                          {FormFactorInfo::Type::F2p, coupl},
                                          {FormFactorInfo::Type::F2n, -coupl},
                                          {FormFactorInfo::Type::FA, coupl}};
        results[{PID::proton(), pid}] = {};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == 23) {
        const std::complex<double> coupl1 = cw*ee*i/(2*sw)-ee*i*sw/(2*cw);
        const std::complex<double> coupl2 = -(cw*ee*i/(2*sw));
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl1},
                                         {FormFactorInfo::Type::F1n, coupl2},
                                         {FormFactorInfo::Type::F2p, coupl1},
                                         {FormFactorInfo::Type::F2n, coupl2},
                                         {FormFactorInfo::Type::FA, coupl2}};
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1n, coupl1},
                                          {FormFactorInfo::Type::F1p, coupl2},
                                          {FormFactorInfo::Type::F2n, coupl1},
                                          {FormFactorInfo::Type::F2p, coupl2},
                                          {FormFactorInfo::Type::FA, coupl2}};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == 22) {
        const std::complex<double> coupl = i*ee;
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                         {FormFactorInfo::Type::F2p, coupl}};
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1n, coupl},
                                          {FormFactorInfo::Type::F2n, coupl}};
        results[{PID::carbon(), pid}] = {{FormFactorInfo::Type::FCoh, 6.0*coupl}};
    } else {
        throw std::runtime_error("LeptonicCurrent: Invalid probe");
    }

    return results;
}

achilles::Currents LeptonicCurrent::CalcCurrents(const std::vector<FourVector> &p,
                                                 const double&) const {
    Currents currents;

    // Setup spinors
    FourVector pU, pUBar;
    if(anti) {
        pUBar = -p[1];
        pU = p.back();
    } else {
        pU = -p[1];
        pUBar = p.back();
    }

    //pU.SetPxPyPzE(0.0,0.0,-620.0,-620.0);
    //pUBar.SetPxPyPzE(-3.68982132e+01, 8.37979251e+00, 1.77567478e+02, 1.81554148e+02);

    spdlog::debug("k = {}, m = {}",pU, pU.M2());
    spdlog::debug("kp = {}, m = {}",pUBar, pUBar.M2());

    std::array<Spinor, 2> ubar, u;
    ubar[0] = UBarSpinor(-1, pUBar);
    ubar[1] = UBarSpinor(1, pUBar);
    u[0] = USpinor(-1, pU);
    u[1] = USpinor(1, pU);

    // Calculate currents
    Current result;
    double q2 = (p[1] - p.back()).M2();
    //double q2 = (pU + pUBar).M2();
    spdlog::debug("Q2 = {}", q2);
    std::complex<double> prop = std::complex<double>(0, 1)/(q2-mass*mass-std::complex<double>(0, 1)*mass*width);
    spdlog::trace("Calculating Current for {}", pid);
    for(size_t i = 0; i < 2; ++i) {
        for(size_t j = 0; j < 2; ++j) {
            std::vector<std::complex<double>> subcur(4);
            for(size_t mu = 0; mu < 4; ++mu) {
                subcur[mu] = ubar[i]*(coupl_left*SpinMatrix::GammaMu(mu)*SpinMatrix::PL()
                                    + coupl_right*SpinMatrix::GammaMu(mu)*SpinMatrix::PR())*u[j]*prop;
                spdlog::trace("Current[{}][{}] = {}", 2*i+j, mu, subcur[mu]);
            }
            result.push_back(subcur);
        }
    }
    currents[pid] = result;

    return currents;
}

void HardScattering::SetProcessGroup(const Process_Group &process_group) {
    spdlog::debug("Adding Process Group: {}", process_group);
    m_group = process_group;
#ifndef ENABLE_BSM
    // TODO: Handle this correctly for different leptonic initial states
    //       This is probably only an issue when including masses,
    //       but those should probably be in a different process group
    m_current.Initialize(process_group.Processes()[0]);
    SMFormFactor = m_current.GetFormFactor();
#endif
}

achilles::Currents HardScattering::LeptonicCurrents(const std::vector<FourVector> &p,
                                                    const double &mu2) const {
#ifdef ENABLE_BSM
    // TODO: Move adapter code into Sherpa interface code
    std::vector<std::array<double, 4>> mom(p.size());
    std::vector<int> pids;
    for(const auto &elm : m_group.Processes()[0].m_mom_map) {
        pids.push_back(static_cast<int>(elm.second));
        mom[elm.first] = (p[elm.first]/1_GeV).Momentum();
        spdlog::debug("PID: {}, Momentum: ({}, {}, {}, {})", pids.back(),
                      mom[elm.first][0], mom[elm.first][1], mom[elm.first][2], mom[elm.first][3]); 
    }
    auto currents = p_sherpa -> Calc(pids, mom, mu2);

    for(auto &current : currents) { 
        spdlog::trace("Current for {}", current.first);
        for(size_t i = 0; i < current.second.size(); ++i) {
            for(size_t j = 0; j < current.second[0].size(); ++j) {
                current.second[i][j] /= pow(1_GeV, static_cast<double>(mom.size())-3);
                spdlog::trace("Current[{}][{}] = {}", i, j, current.second[i][j]);
            }
        }
    }
    return currents;
#else
    return m_current.CalcCurrents(p, mu2);
#endif
}

std::vector<double> HardScattering::CrossSection(Event &event) const {

    //Want to rotate whole system so that q is along z
    auto q = event.Momentum()[1] - event.Momentum().back();
    spdlog::debug("q before Rot = {}", q);
    auto qalongz = q.AlignZ();
    for(size_t i = 0; i < 5; ++i) {
        spdlog::debug("p[{}] before Rot = {}", i, event.Momentum()[i]);
        event.Momentum()[i] = event.Momentum()[i].Rotate(qalongz);
        spdlog::debug("p[{}] after Rot = {}", i, event.Momentum()[i]);
    }
    q = q.Rotate(qalongz);
    spdlog::debug("q after Rot = {}", q);

    // Calculate leptonic currents
    auto leptonCurrent = LeptonicCurrents(event.Momentum(), 100);

    
    std::vector<std::vector<std::complex<double>>> lmunu(4,std::vector<std::complex<double>>(4));
    // Let's compute the leptonic respone tensor for testing purposes
    const size_t nlep_spins = leptonCurrent.begin()->second.size();
    for(size_t i = 0; i  < nlep_spins; ++i) {
        for(size_t mu = 0; mu < 4; ++mu) {
            for(size_t nu = 0; nu < 4; ++nu) {
                for(const auto &lcurrent : leptonCurrent) {
                    lmunu[mu][nu] += std::conj(lcurrent.second[i][mu])*lcurrent.second[i][nu];
                }
            }
        }
    }
    

    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            spdlog::debug("l[{}][{}] = {}",mu,nu,lmunu[mu][nu]);
        }
    }
    


    //std::string filename = "data/pke12_tot.data";

    //static SpectralFunction specf(filename);


    // Calculate the hadronic currents
    // TODO: Clean this up and make generic for the nuclear model
    // TODO: Move this to initialization to remove check each time
    static std::vector<NuclearModel::FFInfoMap> ffInfo;
    if(ffInfo.empty()) {
        ffInfo.resize(m_group.Processes().size());
        for(const auto &current : leptonCurrent) {
            size_t idx = 0;
            for(const auto &process : m_group.Processes()) {
#ifdef ENABLE_BSM
                // TODO: Handle MEC
                ffInfo[idx++][current.first] = p_sherpa -> FormFactors(process.state.first[0], current.first);
#else
                // TODO: Define values somewhere
                ffInfo[idx++][current.first] = SMFormFactor.at({process.state.first[0], current.first});
#endif
            }
        }
    }
    auto hadronCurrent = m_nuclear -> CalcCurrents(event, ffInfo);

    std::vector<double> amps2(hadronCurrent.size());
    //const size_t nlep_spins = leptonCurrent.begin()->second.size();
    const size_t nhad_spins = m_nuclear -> NSpins();
    for(size_t i = 0; i  < nlep_spins; ++i) {
        for(size_t j = 0; j < nhad_spins; ++j) {
            double sign = 1.0;
            std::vector<std::complex<double>> amps(hadronCurrent.size());
            for(size_t mu = 0; mu < 4; ++mu) {
                for(const auto &lcurrent : leptonCurrent) {
                    auto boson = lcurrent.first;
                    for(size_t k = 0; k < hadronCurrent.size(); ++k) {
                        if(hadronCurrent[k].find(boson) != hadronCurrent[k].end()){
                            amps[k] += sign*lcurrent.second[i][mu]*hadronCurrent[k][boson][j][mu];
			             }
                    }
                }
                sign = -1.0;
            }
            for(size_t k = 0; k < hadronCurrent.size(); ++k)
                amps2[k] += std::norm(amps[k]);
        }
    }


    double spin_avg = 1;
    if(!ParticleInfo(m_group.Processes()[0].ids[0]).IsNeutrino()) spin_avg *= 2;
    if(m_nuclear -> NSpins() > 1) spin_avg *= 2;

    static constexpr double to_nb = 1e6;
    std::vector<double> xsecs(m_group.Processes().size());
    for(size_t i = 0; i < m_group.Processes().size(); ++i) {
        // TODO: Handle the case for MEC
        double mass = ParticleInfo(m_group.Processes()[i].state.first[0]).Mass();
        double flux = 2*event.Momentum()[1].E()*2*sqrt(event.Momentum()[0].P2() + mass*mass);

        //double E = -event.Momentum()[0].E() + Constant::mN;
        //double spec = specf(event.Momentum()[0].P(), E);

        FourVector piN_p4 = event.Momentum()[2] + event.Momentum()[3];
        double piN_inv = piN_p4.M();

        if (piN_inv < 1076.957 || piN_inv > 1400.0) amps2[i] = 0.;
        //else amps2[i] = pow(1_MeV,-2);
        spdlog::debug("contraction = {}", amps2[i]);

        //amps2[i] = 1.0;
        // TODO: Correct this flux
        // double flux = 4*sqrt(pow(event.Momentum()[0]*event.Momentum()[1], 2) 
        //                      - event.Momentum()[0].M2()*event.Momentum()[1].M2());


        // Testing phase space only
        
        xsecs[i] = amps2[i]*Constant::HBARC2/spin_avg/flux*to_nb;
        //xsecs[i] = amps2[i]/spin_avg;
        spdlog::debug("Xsec[{}] = {}", i, xsecs[i]);
    }

    for(size_t i = 0; i < 5; ++i) {
        spdlog::debug("p[{}] before Rotback = {}", i, event.Momentum()[i]);
        event.Momentum()[i] = event.Momentum()[i].RotateBack(qalongz);
        spdlog::debug("p[{}] after Rotback = {}", i, event.Momentum()[i]);
    }

    //delete specf;
    return xsecs;
}

bool HardScattering::FillEvent(Event& event, const std::vector<double> &xsecs) const {
    // TODO: Handle multiple different beam particles
    if(!event.CalcTotalCrossSection(xsecs))
        return false;

    auto idx = Random::Instance().SelectIndex(xsecs);
    auto process = m_group.Processes()[idx];
   
    event.InitializeLeptons(process);
    event.InitializeHadrons(process);

    return true;
}
