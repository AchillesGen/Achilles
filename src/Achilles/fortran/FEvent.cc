// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Event.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Particle.hh"

extern "C" {
void GetEventMomentum(achilles::Event *e, size_t i, achilles::FourVector *p) {
    p = &(e->Momentum()[i]);
}

void GetEventNucleus(achilles::Event *e, achilles::Nucleus *nuc) {
    nuc = e->CurrentNucleus().get();
}
}
