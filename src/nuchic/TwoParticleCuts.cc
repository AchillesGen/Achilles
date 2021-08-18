#include "nuchic/TwoParticleCuts.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Units.hh"

bool nuchic::DeltaThetaCut::MakeCut(const FourVector &p1, const FourVector &p2) const {
    return CheckCut(std::acos(p1.CosAngle(p2))/1.0_deg);
}

bool nuchic::InvariantMassCut::MakeCut(const FourVector &p1, const FourVector &p2) const {
    return CheckCut((p1 + p2).M());
}
