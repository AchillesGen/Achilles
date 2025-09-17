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
