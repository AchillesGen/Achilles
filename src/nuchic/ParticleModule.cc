#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include "nuchic/ThreeVector.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Particle.hh"

namespace py = pybind11;

PYBIND11_MODULE(particle, m) {
    py::object vectors = (py::object) py::module::import("vectors");

    py::class_<Particle, std::shared_ptr<Particle>>(m, "Particle")
        // Constructors
        .def(py::init<const int&, const FourVector&, const ThreeVector&,
                      const int&, const std::vector<int>&, const std::vector<int>&>(),
             py::arg("pid") = 0, py::arg("momentum") = FourVector(),
             py::arg("position") = ThreeVector(), py::arg("status") = 0,
             py::arg("mothers") = std::vector<int>(),
             py::arg("daughters") = std::vector<int>())
        // Setters
        .def("set_pid", &Particle::SetPID)
        .def("set_position", &Particle::SetPosition)
        .def("set_momentum", &Particle::SetMomentum)
        .def("set_status", &Particle::SetStatus)
        .def("set_mothers", &Particle::SetMothers)
        .def("set_daughters", &Particle::SetDaughters)
        .def("add_mother", &Particle::AddMother)
        .def("add_daughter", &Particle::AddDaughter)
        .def("set_formation_zone", &Particle::SetFormationZone)
        // Getters
        .def("pid", &Particle::PID)
        .def("position", &Particle::Position)
        .def("momentum", &Particle::Momentum)
        .def("beta", &Particle::Beta)
        .def("status", &Particle::Status)
        .def("mothers", &Particle::Mothers)
        .def("daughters", &Particle::Daughters)
        .def("formation_zone", &Particle::FormationZone)
        .def("mass", &Particle::Mass)
        .def("px", &Particle::Px)
        .def("py", &Particle::Py)
        .def("pz", &Particle::Pz)
        .def("energy", &Particle::E)
        .def("radius", &Particle::Radius)
        // Logical Functions
        .def("is_in_formation_zone", &Particle::InFormationZone)
        .def("is_background", &Particle::IsBackground)
        .def("is_propagating", &Particle::IsPropagating)
        .def("is_final", &Particle::IsFinal)
        // Functions
        .def("propagate", &Particle::Propagate)
        .def("back_propagate", &Particle::BackPropagate)
        // Comparison Operators
        .def(py::self == py::self)
        .def(py::self != py::self)
        // String Methods
        .def("__str__", &Particle::ToString)
        .def("__repr__", &Particle::ToString);

    m.def("check_instances", [](py::list l) {
        return py::make_tuple(
            py::isinstance<Particle>(l[0]),
            py::isinstance<Particle>(l[1])
        );
    });

}
