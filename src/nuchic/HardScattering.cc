#include <iostream>
#include <utility>

#include "nuchic/HardScattering.hh"
#include "nuchic/Constants.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/HardScatteringFactory.hh"
#include "nuchic/Event.hh"
#include "nuchic/Random.hh"

// Aliases for most common types
using nuchic::Particles;
using nuchic::HardScattering;
using nuchic::QESpectral;
using nuchic::QEGlobalFermiGas;

nuchic::HardScattering::HardScattering(RunMode mode) 
    : m_mode{mode} {
    const std::string form_factor_file = "FormFactors.yml";
    YAML::Node node = YAML::LoadFile(form_factor_file);
    spdlog::debug("Setting up form factors");
    const auto vectorFF = node["vector"].as<std::string>();
    const auto axialFF = node["axial"].as<std::string>();
    const auto coherentFF = node["coherent"].as<std::string>();
    m_form_factor = FormFactorBuilder().Vector(vectorFF, node[vectorFF])
                                       .AxialVector(axialFF, node[axialFF])
                                       .Coherent(coherentFF, node[coherentFF])
                                       .build();
}

int HardScattering::LeptonVariables() const {
    switch(m_mode) {
        case RunMode::FixedAngleEnergy:
            return 1;
        case RunMode::FixedAngle:
            return 2;
        case RunMode::FullPhaseSpace:
            return 3;
    }
    return -1;
}

void HardScattering::AddProcess(const nuchic::Process_Info &process) {
    spdlog::debug("Adding Process: {}", process);
    m_leptonicProcesses.push_back(process);
}

nuchic::Currents HardScattering::LeptonicCurrents(const std::vector<FourVector> &p,
                                                  const double &mu2) const {
    std::vector<std::array<double, 4>> mom(p.size());
    auto qVec = p[1];
    for(size_t i = 2; i < p.size()-1; ++i) {
        qVec -= p[i];
    }
    auto rotMat = qVec.AlignZ();
    mom[0] = (-qVec/1_GeV).Rotate(rotMat).Momentum();
    mom[1] = (-p[1]/1_GeV).Rotate(rotMat).Momentum();
    for(size_t i = 2; i < p.size()-1; ++i) {
        mom[i] = (p[i]/1_GeV).Rotate(rotMat).Momentum();
    }
    mom.back() = {};

    size_t idx = 0;
    std::vector<int> pids;
    for(const auto &pid : m_leptonicProcesses[0].m_ids) {
        pids.emplace_back(pid);
        spdlog::debug("PID: {}, Momentum: ({}, {}, {}, {})", pids[idx],
                      mom[idx][0], mom[idx][1], mom[idx][2], mom[idx][3]); 
        ++idx;
    }
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

nuchic::FormFactorArray HardScattering::CouplingsFF(const FormFactor::Values &formFactors,
                                                    const std::vector<FormFactorInfo> &ffInfo) const {
    FormFactorArray results{};

    for(const auto & ff : ffInfo) {
        spdlog::trace("Form Factor: {}, Coupling: {}", ff.form_factor, ff.coupling);
        switch(ff.form_factor) {
            case FormFactorInfo::Type::F1p:
                results[0] += formFactors.F1p*ff.coupling;
                break;
            case FormFactorInfo::Type::F1n:
                results[0] += formFactors.F1n*ff.coupling;
                break;
            case FormFactorInfo::Type::F2p:
                results[1] += formFactors.F2p*ff.coupling;
                break;
            case FormFactorInfo::Type::F2n:
                results[1] += formFactors.F2n*ff.coupling;
                break;
            case FormFactorInfo::Type::FA:
                results[2] += formFactors.FA*ff.coupling;
                break;
            case FormFactorInfo::Type::FCoh:
                results[3] += formFactors.Fcoh*ff.coupling;
                break;
        }
    }

    return results;
}

nuchic::Current HardScattering::CalculateHadronicCurrent(const std::array<Spinor, 2> &ubar,
                                                         const std::array<Spinor, 2> &u,
                                                         const FourVector &qVec,
                                                         const FormFactorArray &ffVal) const {

    Current result;
    std::array<SpinMatrix, 4> gamma{};
    for(size_t mu = 0; mu < 4; ++mu) {
        gamma[mu] = ffVal[0]*SpinMatrix::GammaMu(mu) + ffVal[2]*SpinMatrix::GammaMu(mu)*SpinMatrix::Gamma_5();
        // gamma[mu] = ffVal[2]*SpinMatrix::GammaMu(mu)*SpinMatrix::Gamma_5();
        // gamma[mu] = ffVal[0]*SpinMatrix::GammaMu(mu);
        double sign = 1;
        for(size_t nu = 0; nu < 4; ++nu) {
            gamma[mu] += std::complex<double>(0, 1)*(ffVal[1]*SpinMatrix::SigmaMuNu(mu, nu)*sign*qVec[nu]/(2*Constant::mN));
            sign = -1;
        }
    }

    for(size_t i = 0; i < 2; ++i) {
        for(size_t j = 0; j < 2; ++j) {
            std::vector<std::complex<double>> subcur(4);
            for(size_t mu = 0; mu < 4; ++mu) {
                subcur[mu] = ubar[i]*gamma[mu]*u[j];
            }
            result.push_back(subcur);
        }
    }

    return result;
}

void HardScattering::GeneratePhaseSpace(const std::vector<double> &rans, Event &event) const {
    // Separate the random variables
    std::vector<double> leptonRans(rans.begin(), rans.begin() + LeptonVariables());
    std::vector<double> hadronRans(rans.begin() + LeptonVariables(), rans.end());

    // Generate leptonic phase space point
    GenerateLeptons(leptonRans, event);

    // Generate hadonic phase space point
    FourVector Q = event.Momentum()[0];
    for(const auto &lepton : event.Momentum()) {
        if(lepton == event.Momentum()[0]) continue;
        Q -= lepton;
    }
    GenerateHadrons(hadronRans, Q, event);
}

void HardScattering::GenerateLeptons(const std::vector<double> &leptonRans, Event &event) const {
    double phi = dPhi*leptonRans[0];

    double cosT{}, sinT{}, Elepton{};
    // event.PhaseSpace().weight = dPhi;
    // The following is included in the matrix element
    // leptons[1].E()*leptons[1].Momentum().P();
    switch(m_mode) {
        case RunMode::FixedAngleEnergy:
            Elepton = m_lepton_energy;
            cosT = std::cos(m_angle);
            sinT = std::sin(m_angle);
            break;
        case RunMode::FixedAngle:
            Elepton = event.Momentum()[0].E()*leptonRans[1];
            cosT = std::cos(m_angle);
            sinT = std::sin(m_angle);
            // TEST: fixed energy
            // Elepton = event.PhaseSpace().momentum[0].E()*Constant::mN/(event.PhaseSpace().momentum[0].E()*(1-cosT)+Constant::mN);
            // event.PhaseSpace().weight *= event.PhaseSpace().momentum[0].E();
            break;
        case RunMode::FullPhaseSpace:
            Elepton = event.Momentum()[0].E()*leptonRans[1];
            cosT = dCos*leptonRans[2] - 1;
            sinT = sqrt(1-cosT*cosT);
            // event.PhaseSpace().weight *= dCos*event.PhaseSpace().momentum[0].E();
            break;
    }
    event.Momentum().emplace_back(Elepton,
                                  Elepton*sinT*cos(phi),
                                  Elepton*sinT*sin(phi),
                                  Elepton*cosT);

    // TODO: Move to a better location?
    event.MatrixElements().resize(12);
    for(size_t i = 0; i < event.MatrixElements().size(); ++i) {
        event.MatrixElement(i).inital_state.push_back(PID::electron());
        event.MatrixElement(i).final_state.push_back(PID::electron());
    }
}

bool nuchic::Coherent::InitializeEvent(nuchic::Event &event) {
    // Calculate and store total cross section
    if(!event.TotalCrossSection())
        return false;

    // There is only one matrix element for coherent scattering
    event.InitializeLeptons(0);
    event.InitializeCoherent();

    return true;
}

bool nuchic::Quasielastic::InitializeEvent(nuchic::Event &event) {
    // Calculate and store total cross section
    if(!event.TotalCrossSection())
        return false;

    // Select a matrix element to use
    auto ime = SelectMatrixElement(event);

    // Initialize the leptons for the event
    event.InitializeLeptons(ime);

    // TODO: This is super slow
    // Select the kicked nucleon
    // For the base quasielastic case, we will use global Fermi momentum
    // i.e. we will take MatrixElement(0) is a proton 
    // and MatrixElement(1) is a neutron
    // auto tmp = GetRNG() -> uniform<size_t>(0, event.CurrentNucleus() -> NProtons());
    // size_t idx = ime == 0 ? event.CurrentNucleus() -> ProtonsIDs()[tmp]
    //     : event.CurrentNucleus() -> NeutronsIDs()[tmp];
    // event.InitializeHadrons({{idx, 0, 3}});
    
    event.InitializeHadrons({{ime, 0, Processes()[0].m_ids.size()-1}});

    return true;
}

size_t HardScattering::SelectMatrixElement(nuchic::Event &event) const {
    std::vector<double> probs = event.EventProbs();
    double rand = Random::Instance().Uniform(0.0, 1.0);
    return static_cast<size_t>(std::distance(probs.begin(),
                   std::lower_bound(probs.begin(), probs.end(), rand)))-1;
}

void QESpectral::GenerateHadrons(const std::vector<double> &rans,
                                 const FourVector &Q, Event &event) const {
    // TEST: proton at rest
    // event.PhaseSpace().momentum.emplace_back(Constant::mN, 0, 0, 0);
    // auto pout = event.PhaseSpace().momentum[0] + event.PhaseSpace().momentum[2] - event.PhaseSpace().momentum[1];
    // event.PhaseSpace().momentum.push_back(pout);
    // return;

    static constexpr double dp = 1000; // Hard code the maximum allowed momentum

    // Generate phase space
    double cosT = dCos*rans[0] - 1;
    double sinT = sqrt(1-cosT*cosT);
    double phi = dPhi*rans[1];
    double p = dp*rans[2];

    ThreeVector tmp = { p*sinT*cos(phi), p*sinT*sin(phi), p*cosT };
    tmp += Q.Vec3(); 
    
    // Pauli Blocking
    // TODO: Is this correct? This would assume a global Fermi gas,
    //       but this is then inconsistent with a local Fermi gas for the cascade
    // event.PhaseSpace().weight *= tmp.Magnitude() > 225 ? 1 : 0;

    double Epp = sqrt(pow(Constant::mN, 2)+tmp.P2()); 
    double Ep = Constant::mN+Q.E()-Epp;
    event.Momentum().emplace_back(Ep, p*sinT*cos(phi),
                                  p*sinT*sin(phi),
                                  p*cosT);
    event.Momentum().emplace_back(tmp, Epp);
    // event.PhaseSpace().weight *= Ep > 0 ? dp*p*p*dCos*dPhi : 0;
}

void QEGlobalFermiGas::GenerateHadrons(const std::vector<double> &rans,
                                       const FourVector &Q, Event &event) const {
    static const double dp = 1000; // Hard code the maximum allowed momentum
    static constexpr double Ep = 20.0;

    // Generate phase space
    double phi = 2*M_PI*rans[0];
    double p = dp*rans[1];
    double Ef = sqrt(pow(p, 2)+pow(Constant::mN, 2));
    
    double arg1 = pow(Q.E()-Ep+Ef,2);
    double arg2 = pow(p, 2)+pow(Constant::mN, 2)+pow(Q.P(), 2);
	    
    double cosT = (arg1-arg2)/(2*p*Q.P());
    double sinT = sqrt(1-cosT*cosT);
    if(std::abs(cosT) > 1) {
        // event.PhaseSpace().weight = 0;
        return;
    } 

    ThreeVector tmp = { p*sinT*cos(phi), p*sinT*sin(phi), p*cosT };
    tmp += Q.Vec3();
    
    // Pauli Blocking
    
    double Epp = sqrt(pow(Constant::mN, 2)+tmp.P2()); 
    
    event.Momentum().emplace_back(Ep, p*sinT*cos(phi),
                                             p*sinT*sin(phi),
                                             p*cosT);
    event.Momentum().emplace_back(tmp, Epp);
    // event.PhaseSpace().weight *= dp*p*p*dPhi;
}
