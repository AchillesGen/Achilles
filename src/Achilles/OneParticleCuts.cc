#include "Achilles/OneParticleCuts.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Units.hh"

bool achilles::EnergyCut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.E());
}

bool achilles::MomentumCut::MakeCut(const FourVector &mom) const {
    spdlog::trace("Momentum = {}, Cut Result = {}", mom.P(), CheckCut(mom.P()));
    return CheckCut(mom.P());
}

bool achilles::AngleThetaCut::MakeCut(const FourVector &mom) const {
    spdlog::trace("Theta = {}, Cut Result = {}", mom.Theta() / 1.0_deg,
                  CheckCut(mom.Theta() / 1.0_deg));
    return CheckCut(mom.Theta() / 1.0_deg);
}

bool achilles::TransverseMomentumCut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.Pt());
}

bool achilles::ETheta2Cut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.E() * pow(mom.Theta(), 2));
}

achilles::Q2Cut::Q2Cut(const YAML::Node &node)
    : OneParticleCut(node), m_pid{node["BeamId"].as<PID>()}, m_beam{node["Beam"].as<double>()} {}

bool achilles::Q2Cut::MakeCut(const FourVector &mom) const {
    double mass = ParticleInfo(m_pid).Mass();
    FourVector beam{m_beam, 0, 0, sqrt(m_beam * m_beam - mass * mass)};
    double Q2 = -(beam - mom).M2();
    return CheckCut(Q2);
}
