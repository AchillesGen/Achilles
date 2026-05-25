// SPDX-FileCopyrightText: 2018-2024 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/PyBindings.hh"

PYBIND11_MODULE(_achilles, m) {
    // Utilities
    py::module utilities = m.def_submodule("utilities", "achilles utilities");
    LoggingModule(utilities);
    ConstantsModule(utilities);
    InterpolationModule(utilities);

    // Physics Objects
    py::module physics = m.def_submodule("physics", "achilles physics objects");
    VectorModule(physics);
    ParticleInfoModule(physics);
    ParticleModule(physics);
    NucleusModule(physics);

    // Calculation Objects
    InteractionsModule(m);
    CascadeModule(m);
}
