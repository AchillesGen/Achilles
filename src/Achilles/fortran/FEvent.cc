#include "Achilles/Event.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Particle.hh"

extern "C" {
    void GetEventMomentum(achilles::Event *e, size_t i, achilles::FourVector *p) {
        p = &(e -> Momentum()[i]);
    }

    void SetEventMatrixWgt(achilles::Event *e, size_t i, double wgt) {
        e -> MatrixElementWgt(i) = wgt;
    }

    bool EventTotalCrossSection(achilles::Event *e) {
        return e -> TotalCrossSection();
    }

    size_t EventSelectNucleon(achilles::Event *e) {
        return e -> SelectNucleon();
    }

    void GetEventNucleus(achilles::Event *e, achilles::Nucleus *nuc) {
        nuc = e -> CurrentNucleus().get();
    }

    void SetEventMEWgt(achilles::Event *e, double wgt) {
        e -> SetMEWeight(wgt);
    }
}
