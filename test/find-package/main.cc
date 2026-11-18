// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/FourVector.hh"
#include "Achilles/Interactions.hh"

#include <iostream>

int main() {
    achilles::FourVector p{1, 0, 0, 1};
    if(p.M2() != 0) return 1;

    // Cascade models self-register at static init, so this fails if the linker drops them.
    if(!achilles::InteractionFactory::IsRegistered("ConstantInteraction")) {
        std::cerr << "Cascade interactions were not registered\n";
        return 1;
    }
    return 0;
}
