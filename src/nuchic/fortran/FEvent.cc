#include "nuchic/Event.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Particle.hh"

extern "C" {
    void GetEventMomentum(nuchic::Event *e, size_t i, nuchic::FourVector *p) {
        p = &(e -> Momentum()[i]);
    }

    void SetEventMatrixWgt(nuchic::Event *e, size_t i, double wgt) {
        e -> MatrixElementWgt(i) = wgt;
    }

    bool EventTotalCrossSection(nuchic::Event *e) {
        return e -> TotalCrossSection();
    }

    size_t EventSelectNucleon(nuchic::Event *e) {
        return e -> SelectNucleon();
    }

    void GetEventNucleus(nuchic::Event *e, nuchic::Nucleus *nuc) {
        nuc = e -> CurrentNucleus().get();
    }

    void SetEventMEWgt(nuchic::Event *e, double wgt) {
        e -> SetMEWeight(wgt);
    }
}
