#include "Achilles/Event.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Particle.hh"

extern "C" {
void GetEventMomentum(achilles::Event *e, size_t i, achilles::FourVector *p) {
    p = &(e->Momentum()[i]);
}
}
