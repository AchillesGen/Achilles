#include "nuchic/OneParticleCuts.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Units.hh"

bool nuchic::EnergyCut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.E());
}

bool nuchic::MomentumCut::MakeCut(const FourVector &mom) const {
    spdlog::trace("Momentum = {}, Cut Result = {}", mom.P(), CheckCut(mom.P()));
    return CheckCut(mom.P());
}

bool nuchic::AngleThetaCut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.Theta()/1.0_deg);
}

bool nuchic::TransverseMomentumCut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.Pt());
}

bool nuchic::ETheta2Cut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.E()*pow(mom.Theta(), 2));
}
