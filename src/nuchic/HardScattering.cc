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
    : m_mode{mode} {}

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

std::vector<double> HardScattering::LeptonicTensor(const std::vector<FourVector> &p,
                                                   const double &mu2) const {
    std::vector<std::array<double, 4>> mom;
    for(const auto &pi : p) mom.emplace_back((pi/1_GeV).Momentum());
    std::swap(mom[1], mom[2]);
    return p_sherpa -> Calc(m_pids, mom, mu2);
}

void HardScattering::GeneratePhaseSpace(const std::vector<double> &rans, Event &event) const {
    // Separate the random variables
    std::vector<double> leptonRans(rans.begin(), rans.begin() + LeptonVariables());
    std::vector<double> hadronRans(rans.begin() + LeptonVariables(), rans.end());

    // Generate leptonic phase space point
    GenerateLeptons(leptonRans, event);

    // Generate hadonic phase space point
    FourVector Q = event.PhaseSpace().momentum[0];
    for(const auto &lepton : event.PhaseSpace().momentum) {
        if(lepton == event.PhaseSpace().momentum[0]) continue;
        Q -= lepton;
    }
    GenerateHadrons(hadronRans, Q, event);
}

void HardScattering::GenerateLeptons(const std::vector<double> &leptonRans, Event &event) const {
    double phi = dPhi*leptonRans[0];

    double cosT{}, sinT{}, Elepton{};
    event.PhaseSpace().weight = dPhi;
    // The following is included in the matrix element
    // leptons[1].E()*leptons[1].Momentum().P();
    switch(m_mode) {
        case RunMode::FixedAngleEnergy:
            Elepton = m_lepton_energy;
            cosT = std::cos(m_angle);
            sinT = std::sin(m_angle);
            break;
        case RunMode::FixedAngle:
            Elepton = event.PhaseSpace().momentum[0].E()*leptonRans[1];
            cosT = std::cos(m_angle);
            sinT = std::sin(m_angle);
            event.PhaseSpace().weight *= event.PhaseSpace().momentum[0].E();
            break;
        case RunMode::FullPhaseSpace:
            Elepton = event.PhaseSpace().momentum[0].E()*leptonRans[1];
            cosT = dCos*leptonRans[2] - 1;
            sinT = sqrt(1-cosT*cosT);
            event.PhaseSpace().weight *= dCos*event.PhaseSpace().momentum[0].E();
            break;
    }
    event.PhaseSpace().momentum.emplace_back(Elepton*sinT*cos(phi),
                                             Elepton*sinT*sin(phi),
                                             Elepton*cosT, Elepton);

    // TODO: Move to a better location?
    event.MatrixElements().resize(12);
    for(size_t i = 0; i < event.MatrixElements().size(); ++i) {
        event.MatrixElement(i).inital_state.push_back(PID::electron());
        event.MatrixElement(i).final_state.push_back(PID::electron());
    }
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
    // event.InitializeHadrons({{idx, 2, 3}});
    
    event.InitializeHadrons({{ime, 2, 3}});

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
    static const double dp = 1000; // Hard code the maximum allowed momentum

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
    event.PhaseSpace().momentum.emplace_back(p*sinT*cos(phi),
                                             p*sinT*sin(phi),
                                             p*cosT, Ep);
    event.PhaseSpace().momentum.emplace_back(tmp, Epp);
    event.PhaseSpace().weight *= Ep > 0 ? dp*p*p*dCos*dPhi : 0;
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
        event.PhaseSpace().weight = 0;
        return;
    } 

    ThreeVector tmp = { p*sinT*cos(phi), p*sinT*sin(phi), p*cosT };
    tmp += Q.Vec3();
    
    // Pauli Blocking
    
    double Epp = sqrt(pow(Constant::mN, 2)+tmp.P2()); 
    
    event.PhaseSpace().momentum.emplace_back(p*sinT*cos(phi),
                                             p*sinT*sin(phi),
                                             p*cosT, Ep);
    event.PhaseSpace().momentum.emplace_back(tmp, Epp);
    event.PhaseSpace().weight *= dp*p*p*dPhi;
}
