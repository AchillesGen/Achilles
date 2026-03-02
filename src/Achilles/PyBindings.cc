#include "Achilles/PyBindings.hh"
#include "Achilles/AchillesInit.hh"
#include "Achilles/Version.hh"
#include "Achilles/fortran/FNuclearModel.hh"
#include "fmt/core.h"

PYBIND11_MODULE(_achilles, m) {
    // Achilles Binding
    std::vector<std::string> shargs;
    m.attr("__version__") = py::cast(fmt::format("achilles {}", ACHILLES_VERSION));
    m.def("set_paths", &achilles::InitializePaths, "Initializes the search paths for Achilles.");
    m.def("register_fortran_models", &achilles::FortranModel::RegisterModels,
          "Register all the fortran models");
    m.def("initialize_logging", &achilles::InitializeLogging,
          "Sets logger settings. Can only be called once, before generating events.",
          py::arg("verbosity") = 2, py::arg("log_verbosity") = 2,
          py::arg("logfile") = "achilles.log");
    m.def("generate_events", &achilles::GenerateEvents, "Runs the event generator.",
          py::arg("runcard") = "run.yml", py::arg("sherpa_args") = shargs,
          py::arg("batch_mode") = false);
    // verbosity levels:
    // trace = 0
    // debug = 1
    // info = 2
    // warn = 3
    // error = 4
	py::module verbosity=m.def_submodule("verbosity","Logging Verbosity Levels");
	verbosity.attr("trace")=0;
	verbosity.attr("debug")=1;
	verbosity.attr("info")=2;
	verbosity.attr("warn")=3;
	verbosity.attr("error")=4;

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
