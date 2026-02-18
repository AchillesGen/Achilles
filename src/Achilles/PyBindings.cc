#include "Achilles/PyBindings.hh"
#include "Achilles/AchillesInit.hh"
#include "Achilles/Version.hh"
#include "fmt/core.h"

PYBIND11_MODULE(_achilles, m) {
    // Achilles Binding
    std::vector<std::string> shargs;
    m.attr("__version__") = py::cast(fmt::format("achilles {}", ACHILLES_VERSION));
    m.def("set_paths", &achilles::InitializePaths, "Initializes the search paths for Achilles.");
    m.def("generate_events", &achilles::GenerateEvents, "Runs the event generator.",
          py::arg("runcard") = "run.yml", py::arg("shargs") = shargs, py::arg("verbosity") = 0,
          py::arg("log_verbosity") = 0, py::arg("logfile") = "achilles.log",
          py::arg("batch_mode") = false);

    // TODO: This is an outdated interface. We should update or remove
    // // Utilities
    // py::module utilities = m.def_submodule("utilities", "achilles utilities");
    // LoggingModule(utilities);
    // ConstantsModule(utilities);
    // InterpolationModule(utilities);

    // // Physics Objects
    // py::module physics = m.def_submodule("physics", "achilles physics objects");
    // VectorModule(physics);
    // ParticleInfoModule(physics);
    // ParticleModule(physics);
    // NucleusModule(physics);

    // // Calculation Objects
    // InteractionsModule(m);
    // CascadeModule(m);
}
