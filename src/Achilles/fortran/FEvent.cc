#include "Achilles/Event.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Particle.hh"

extern "C" {
void GetEventMomentum(achilles::Event *e, size_t i, achilles::FourVector *p) {
	if(i<e->LeptonsIn().size()) {
		p=&e->LeptonsIn()[i].Momentum();
		return;
	}
	i-=e->LeptonsIn().size();
	if(i<e->HadronsIn().size()) {
		p=&e->HadronsIn()[i].Momentum();
		return;
	}
	i-=e->HadronsIn().size();
	if(i<e->LeptonsOut().size()) {
		p=&e->LeptonsOut()[i].Momentum();
		return;
	}
	i-=e->LeptonsOut().size();
	i-=e->LeptonsIn().size();
	if(i<e->HadronsOut().size()) {
		p=&e->HadronsOut()[i].Momentum();
		return;
	}
	i-=e->HadronsOut().size();
	p=&e->Spectators()[i].Momentum();
}

void GetEventNucleus(achilles::Event *e, achilles::Nucleus *nuc) {
    nuc = e->CurrentNucleus().get();
}
}
