#include "nuchic/HardScattering.hh"
#include "nuchic/Constants.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"

// Aliases for most common types
using nuchic::Particles;
using nuchic::HardScattering;
using nuchic::QESpectral;

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

    std::vector<double> beamRans(rans.begin(), rans.begin() + GetBeam().NVariables());
    std::vector<double> leptonRans(rans.begin()+GetBeam().NVariables(), rans.end());
    leptons[0].SetMomentum(GetBeam().Flux(leptons[0].ID(), beamRans));

    double cosT = 2*leptonRans[0] - 1;
    double sinT = sqrt(1-cosT*cosT);
    double phi = 2*M_PI*leptonRans[1];
    double Elepton = leptons[0].E()*leptonRans[2];
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
    return leptons[0].E()*leptons[1].E()*leptons[1].Momentum().P()*4*M_PI;
}

Particles QESpectral::GenerateHadrons(const std::vector<double> &rans, const FourVector &Q) const {
    Particles hadrons = { Particle(PID::neutron()), Particle(PID::neutron()) };

    double cosT = 2*rans[0] - 1;
    double sinT = sqrt(1-cosT*cosT);
    double phi = 2*M_PI*rans[1];
    double p = GetNucleus() -> FermiMomentum(0)*rans[2];

    ThreeVector tmp = { p*sinT*cos(phi), p*sinT*sin(phi), p*cosT };
    tmp += Q.Vec3();

    double Epp = sqrt(pow(hadrons[1].Mass(), 2)+tmp.P2()); 

    double Ep = hadrons[0].Mass()+Q.E()-Epp;
    Ep = Ep > 0 ? Ep : 0;

    hadrons[0].SetMomentum({p*sinT*cos(phi),
                            p*sinT*sin(phi),
                            p*cosT, Ep});
    hadrons[1].SetMomentum({tmp, Epp});

    return hadrons;
}

double QESpectral::HadronWeight(const Particles &hadrons) const {
    return GetNucleus() -> FermiMomentum(0)*hadrons[0].Momentum().P()*hadrons[0].E()*4*M_PI;
}

double HardScattering::Test(const std::vector<double> &rans, const double &wgt) {
    auto particles = GeneratePhaseSpace(rans);
    auto weight = PhaseSpaceWeight(particles);
    if(m_fill) hist.Fill(particles[2].E(), weight*wgt);

    return weight;
}
