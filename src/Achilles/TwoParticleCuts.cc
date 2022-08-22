#include "Achilles/TwoParticleCuts.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Units.hh"

bool achilles::DeltaThetaCut::MakeCut(const FourVector &p1, const FourVector &p2) const {
    return CheckCut(std::acos(p1.CosAngle(p2))/1.0_deg);
}

bool achilles::InvariantMassCut::MakeCut(const FourVector &p1, const FourVector &p2) const {
    return CheckCut((p1 + p2).M());
}
