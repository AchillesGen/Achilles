#include "Achilles/Achilles.hh"
#include "Achilles/PyBindings.hh"
#include "Achilles/Version.hh"
#include "fmt/core.h"

PYBIND11_MODULE(_achilles, m) {

	// Achilles Binding
	std::vector<std::string> shargs;
	m.attr("version")=py::cast(fmt::format("achilles {}",ACHILLES_VERSION));
	m.def("generate_events",&achilles::GenerateEvents,"Runs the event generator.",
			py::arg("runcard")="run.yml",
			py::arg("shargs")=shargs,
			py::arg("verbosity")=0,
			py::arg("log_verbosity")=0,
			py::arg("logfile")="achilles.log",
			py::arg("batch_mode")=false);

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
