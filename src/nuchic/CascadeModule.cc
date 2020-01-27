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

PYBIND11_MODULE(cascade, m) {
    py::object vectors = (py::object) py::module::import("vectors");
    py::object particle = (py::object) py::module::import("particle");
    py::object nucleus = (py::object) py::module::import("nucleus");
    py::object interactions = (py::object) py::module::import("interactions");

    py::class_<Cascade, std::shared_ptr<Cascade>>(m, "Cascade")
        // Constructors
        .def(py::init<const std::shared_ptr<Interactions>, const double&>(),
                py::arg("interactions"), py::arg("distance") = 0.05)
        // Functions
        .def("kick", &Cascade::Kick)
        .def("reset", &Cascade::Reset)
        .def("set_kicked", &Cascade::SetKicked)
        .def("call", &Cascade::operator())
        .def("__call__", &Cascade::operator(),
                py::arg("particles"), py::arg("fermi_momentum"),
                py::arg("radius2"), py::arg("max_steps") = 10000);
}
