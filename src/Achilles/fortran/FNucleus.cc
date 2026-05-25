// SPDX-FileCopyrightText: 2018-2024 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Nucleus.hh"

extern "C" {
void GetNucleon(achilles::Nucleus *nuc, size_t i, achilles::Particle *part) {
    part = &(nuc->Nucleons()[i]);
}

void AddParticle(achilles::Nucleus *nuc, achilles::Particle *part) {
    achilles::Particle tmp(*part);
    nuc->Nucleons().push_back(tmp);
}
}
