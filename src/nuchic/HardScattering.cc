#include "nuchic/HardScattering.hh"
#include "nuchic/Constants.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/HardScatteringFactory.hh"
#include <iostream>

// Aliases for most common types
using nuchic::Particles;
using nuchic::HardScattering;
using nuchic::QESpectral;
using nuchic::QEGlobalFermiGas;

nuchic::HardScattering::HardScattering(std::shared_ptr<Beam> leptonBeam,
                                       std::shared_ptr<Nucleus> nuc,
                                       RunMode mode)
        : m_leptonBeam{std::move(leptonBeam)}, m_nuc{std::move(nuc)}, m_mode{mode} {
    m_rng = std::make_shared<randutils::mt19937_rng>();
}

int HardScattering::LeptonVariables() const {
    int nvars = GetBeam() -> NVariables();
    switch(m_mode) {
        case RunMode::FixedAngle:
            nvars += 2;
            break;
        case RunMode::FullPhaseSpace:
            nvars += 3;
            break;
    }
    return nvars;
}

double HardScattering::GeneratePhaseSpace(Particles &parts, const std::vector<double> &rans) const {
    std::vector<double> leptonRans(rans.begin(), rans.begin() + LeptonVariables());
    std::vector<double> hadronRans(rans.begin() + LeptonVariables(), rans.end());
    Particles leptons{};
    auto leptonWgt = GenerateLeptons(leptons, leptonRans);

    FourVector Q = leptons[0].Momentum();
    for(const auto &lepton : leptons) {
        if(lepton == leptons[0]) continue;
        Q -= lepton.Momentum();
    }

    auto hadronWgt = GenerateHadrons(parts, hadronRans, Q);
    parts.insert(parts.end(), leptons.begin(), leptons.end());

    return leptonWgt*hadronWgt;
}

double HardScattering::GenerateLeptons(Particles &leptons, const std::vector<double> &rans) const {
    leptons = { Particle(PID::electron()), Particle(PID::electron()) };

    std::vector<double> beamRans(rans.begin(), rans.begin() + GetBeam()->NVariables());
    std::vector<double> leptonRans(rans.begin()+GetBeam()->NVariables(), rans.end());
    leptons[0].SetMomentum(GetBeam()->Flux(leptons[0].ID(), beamRans));

    double phi = dPhi*leptonRans[0];
    double Elepton = leptons[0].E()*leptonRans[1];

    double cosT{}, sinT{};
    double wgt = leptons[0].E()*dPhi;//leptons[1].E()*leptons[1].Momentum().P();
    switch(m_mode) {
        case RunMode::FixedAngle:
            cosT = std::cos(m_angle);
            sinT = std::sin(m_angle); 
            break;
        case RunMode::FullPhaseSpace:
            cosT = dCos*leptonRans[2] - 1;
            sinT = sqrt(1-cosT*cosT);
            wgt *= dCos;
            break;
    }
    leptons[1].SetMomentum({Elepton*sinT*cos(phi),
                            Elepton*sinT*sin(phi),
                            Elepton*cosT, Elepton});
    leptons[0].SetStatus(ParticleStatus::initial_state);
    leptons[1].SetStatus(ParticleStatus::final_state);

    return wgt;
}

double QESpectral::GenerateHadrons(Particles& hadrons, const std::vector<double> &rans, const FourVector &Q) const {
    static const double dp = 4*GetNucleus() -> FermiMomentum(0); 

    // Setup initial state
    auto pid = 2*rans.back() < 1.0 ? PID::proton() : PID::neutron();
    std::vector<std::size_t> indices;
    for(std::size_t i = 0; i < hadrons.size(); ++i) {
        if(hadrons[i].ID() == pid) indices.push_back(i);
    }
    auto initialIdx = GetRNG() -> pick(indices);
    auto initialState = &hadrons[initialIdx];
    initialState -> SetStatus(ParticleStatus::initial_state);
    hadrons.emplace_back(*initialState);
    auto finalState = &hadrons.back();
    finalState -> SetStatus(ParticleStatus::propagating);

    // Generate phase space
    double cosT = dCos*rans[1] - 1;
    double sinT = sqrt(1-cosT*cosT);
    double phi = dPhi*rans[2];
    double p = dp*rans[3];

    ThreeVector tmp = { p*sinT*cos(phi), p*sinT*sin(phi), p*cosT };
    tmp += Q.Vec3();

    double Epp = sqrt(pow(finalState -> Mass(), 2)+tmp.P2()); 
    double Ep = hadrons[initialIdx].Mass()+Q.E()-Epp;
    hadrons[initialIdx].SetMomentum({p*sinT*cos(phi),
                                     p*sinT*sin(phi),
                                     p*cosT, Ep});
    finalState -> SetMomentum({tmp, Epp});
    return dp*pow(hadrons[initialIdx].Momentum().P(), 2)*dCos*dPhi*nNucleonTypes;
}

double QEGlobalFermiGas::GenerateHadrons(Particles &hadrons, const std::vector<double> &rans, const FourVector &Q) const {
    static const double dp = GetNucleus() -> FermiMomentum(0);
    static constexpr double Ep = 20.0;

    // Setup initial state
    auto pid = 2*rans.back() < 1.0 ? PID::proton() : PID::neutron();
    std::vector<std::size_t> indices;
    for(std::size_t i = 0; i < hadrons.size(); ++i) {
        if(hadrons[i].ID() == pid) indices.push_back(i);
    }
    auto initialIdx = GetRNG() -> pick(indices);
    auto initialState = &hadrons[initialIdx];
    initialState -> SetStatus(ParticleStatus::initial_state);
    hadrons.emplace_back(*initialState);
    auto finalState = &hadrons.back();
    finalState -> SetStatus(ParticleStatus::propagating);

    // Generate phase space
    double phi = 2*M_PI*rans[1];
    double p = dp*rans[2];
    double Ef=sqrt(pow(p,2)+pow(initialState -> Mass(),2));
    
    double arg1=pow(Q.E()-Ep+Ef,2);
    double arg2=pow(p,2)+pow(finalState -> Mass(),2)+pow(Q.P(),2);
	    
    double cosT=(arg1-arg2)/(2*p*Q.P());
    if(std::abs(cosT) > 1.0) return 0;
    double sinT = sqrt(1-cosT*cosT);

    ThreeVector tmp = { p*sinT*cos(phi), p*sinT*sin(phi), p*cosT };
    tmp += Q.Vec3();
    double Epp = sqrt(pow(hadrons[1].Mass(), 2)+tmp.P2()); 

    initialState -> SetMomentum({p*sinT*cos(phi),
                                 p*sinT*sin(phi),
                                 p*cosT, Ep});
    finalState -> SetMomentum({tmp, Epp});

    return dp*pow(initialState -> Momentum().P(), 2)*dCos*dPhi*nNucleonTypes;
}
