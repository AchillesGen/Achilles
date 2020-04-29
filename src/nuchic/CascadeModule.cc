#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"
#include "pybind11/functional.h"

#include "nuchic/Cascade.hh"
#include "nuchic/Interactions.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"

namespace py = pybind11;

// These are for convenience
using nuchic::Cascade; 
using nuchic::Interactions;

PYBIND11_MODULE(cascade, m) {
    py::object vectors = static_cast<py::object>(py::module::import("vectors"));
    py::object particle = static_cast<py::object>(py::module::import("particle"));
    py::object nucleus = static_cast<py::object>(py::module::import("nucleus"));
    py::object interactions = static_cast<py::object>(py::module::import("interactions"));

    py::class_<Cascade, std::shared_ptr<Cascade>> cascade(m, "Cascade");

        // Constructors
        cascade.def(py::init<const std::shared_ptr<Interactions>, const Cascade::ProbabilityType&, 
                    const double&>(),
                    py::arg("interactions"), py::arg("prob"), py::arg("distance") = 0.05)
        // Functions
        .def("kick", &Cascade::Kick)
        .def("reset", &Cascade::Reset)
        .def("set_kicked", &Cascade::SetKicked)
        .def("evolve", &Cascade::Evolve,
             py::arg("nucleus"), py::arg("max_steps") = 10000)
	.def("nuwro", &Cascade::NuWro,
             py::arg("nucleus"), py::arg("max_steps") = 10000)
        .def("mean_free_path", &Cascade::MeanFreePath,
             py::arg("nucleus"), py::arg("max_steps") = 10000)
	 .def("mean_free_path_nuwro", &Cascade::MeanFreePath_NuWro,
             py::arg("nucleus"), py::arg("max_steps") = 10000);

    py::enum_<Cascade::ProbabilityType>(cascade, "Probability")
        .value("Gaussian", Cascade::ProbabilityType::Gaussian)
        .value("Pion", Cascade::ProbabilityType::Pion)
        .export_values();
}
