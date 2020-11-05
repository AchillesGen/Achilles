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
    dp = 4*m_nuc -> FermiMomentum(0);
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

Particles HardScattering::GeneratePhaseSpace(const std::vector<double> &rans) const {
    std::vector<double> leptonRans(rans.begin(), rans.begin() + LeptonVariables());
    std::vector<double> hadronRans(rans.begin() + LeptonVariables(), rans.end());
    auto leptons = GenerateLeptons(leptonRans);

    FourVector Q = leptons[0].Momentum();
    for(const auto &lepton : leptons) {
        if(lepton == leptons[0]) continue;
        Q -= lepton.Momentum();
    }

    auto hadrons = GenerateHadrons(hadronRans, Q);

    Particles particles;
    particles.reserve(leptons.size() + hadrons.size());
    particles.insert(particles.end(), leptons.begin(), leptons.end());
    particles.insert(particles.end(), hadrons.begin(), hadrons.end());

    return particles;
}

Particles HardScattering::GenerateLeptons(const std::vector<double> &rans) const {
    Particles leptons = { Particle(PID::electron()), Particle(PID::electron()) };

    std::vector<double> beamRans(rans.begin(), rans.begin() + GetBeam()->NVariables());
    std::vector<double> leptonRans(rans.begin()+GetBeam()->NVariables(), rans.end());
    leptons[0].SetMomentum(GetBeam()->Flux(leptons[0].ID(), beamRans));

    double phi = dPhi*leptonRans[0];
    double Elepton = leptons[0].E()*leptonRans[1];

    double cosT{}, sinT{};
    switch(m_mode) {
        case RunMode::FixedAngle:
            cosT = std::cos(m_angle);
            sinT = std::sin(m_angle); 
            break;
        case RunMode::FullPhaseSpace:
            cosT = dCos*leptonRans[2] - 1;
            sinT = sqrt(1-cosT*cosT);
            break;
    }
    leptons[1].SetMomentum({Elepton*sinT*cos(phi),
                            Elepton*sinT*sin(phi),
                            Elepton*cosT, Elepton});

    return leptons;
}

double HardScattering::PhaseSpaceWeight(const Particles &particles) const {
    Particles leptons(particles.begin(), particles.begin() + 2);
    Particles hadrons(particles.begin()+2, particles.end());

    return LeptonWeight(leptons) * HadronWeight(hadrons);
}

double HardScattering::LeptonWeight(const Particles &leptons) const {
    double wgt = leptons[0].E()*dPhi;//leptons[1].E()*leptons[1].Momentum().P();
    switch(m_mode) {
        case RunMode::FixedAngle:
            break;
        case RunMode::FullPhaseSpace:
            wgt *= dCos;
            break;
    }
    return wgt;
}

Particles QESpectral::GenerateHadrons(const std::vector<double> &rans, const FourVector &Q) const {
    auto pid = 2*rans[0] < 1.0 ? PID::proton() : PID::neutron();
    Particles hadrons = { Particle(pid), Particle(pid) };

    double cosT = dCos*rans[1] - 1;
    double sinT = sqrt(1-cosT*cosT);
    double phi = dPhi*rans[2];
    double p = dp*rans[3];

    ThreeVector tmp = { p*sinT*cos(phi), p*sinT*sin(phi), p*cosT };
    tmp += Q.Vec3();

    double Epp = sqrt(pow(hadrons[1].Mass(), 2)+tmp.P2()); 
    double Ep = hadrons[0].Mass()+Q.E()-Epp;
    hadrons[0].SetMomentum({p*sinT*cos(phi),
                            p*sinT*sin(phi),
                            p*cosT, Ep});
    hadrons[1].SetMomentum({tmp, Epp});

    return hadrons;
}

double QESpectral::HadronWeight(const Particles &hadrons) const {
    return dp*pow(hadrons[0].Momentum().P(), 2)*dCos*dPhi*nNucleonTypes;
}

Particles QEGlobalFermiGas::GenerateHadrons(const std::vector<double> &rans, const FourVector &Q) const {
    auto pid = 2*rans[0] < 1.0 ? PID::proton() : PID::neutron();

    Particles hadrons = { Particle(pid), Particle(pid) };

    double phi = 2*M_PI*rans[1];
    double p = dp*rans[2];
    double Ep = 20.0;
    double Ef=sqrt(pow(p,2)+pow(hadrons[0].Mass(),2));
    
    double arg1=pow(Q.E()-Ep+Ef,2);
    double arg2=pow(p,2)+pow(hadrons[0].Mass(),2)+pow(Q.P(),2);
	    
    double cosT=(arg1-arg2)/(2*p*Q.P());
    if(std::abs(cosT) > 1.0) Ep=-10;
    double sinT = sqrt(1-cosT*cosT);

    ThreeVector tmp = { p*sinT*cos(phi), p*sinT*sin(phi), p*cosT };
    tmp += Q.Vec3();
    double Epp = sqrt(pow(hadrons[1].Mass(), 2)+tmp.P2()); 

    hadrons[0].SetMomentum({p*sinT*cos(phi),
                            p*sinT*sin(phi),
                            p*cosT, Ep});
    hadrons[1].SetMomentum({tmp, Epp});

    return hadrons;
}

double QEGlobalFermiGas::HadronWeight(const Particles &hadrons) const {
    return dp*pow(hadrons[0].Momentum().P(), 2)*dCos*dPhi*nNucleonTypes;
}

double HardScattering::Test(const std::vector<double> &rans, const double &wgt) {
    auto particles = GeneratePhaseSpace(rans);
    auto weight = PhaseSpaceWeight(particles);
    if(m_fill) hist.Fill(particles[2].E(), weight*wgt);

    return weight;
}
