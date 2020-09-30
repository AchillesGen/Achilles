#include "nuchic/PyBindings.hh"

PYBIND11_MODULE(_nuchic, m) {
    // Utilities
    py::module utilities = m.def_submodule("utilities", "nuchic utilities");
    LoggingModule(utilities);
    ConstantsModule(utilities);
    InterpolationModule(utilities);

    // Physics Objects
    py::module physics = m.def_submodule("physics", "nuchic physics objects");
    VectorModule(physics);
    ParticleInfoModule(physics);
    ParticleModule(physics);
    NucleusModule(physics);

    // Calculation Objects
    InteractionsModule(m);
    CascadeModule(m);
}
