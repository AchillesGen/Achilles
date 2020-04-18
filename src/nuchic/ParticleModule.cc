#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include "nuchic/ThreeVector.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Particle.hh"

namespace py = pybind11;

// These are for convenience
using nuchic::Particle;
using nuchic::FourVector;
using nuchic::ThreeVector;

PYBIND11_MODULE(particle, m) {
    py::object vectors = static_cast<py::object>(py::module::import("vectors"));

    py::enum_<nuchic::ParticleStatus>(m, "ParticleStatus")
        .value("internal_test", nuchic::ParticleStatus::internal_test)
        .value("external_test", nuchic::ParticleStatus::external_test)
        .value("propagating", nuchic::ParticleStatus::propagating)
        .value("background", nuchic::ParticleStatus::background)
        .value("escaped", nuchic::ParticleStatus::escaped);

    py::class_<Particle, std::shared_ptr<Particle>>(m, "Particle")
        // Constructors
        .def(py::init<const int&, const FourVector&, const ThreeVector&,
                      const int&, const std::vector<int>&, const std::vector<int>&>(),
             py::arg("pid") = 0, py::arg("momentum") = FourVector(),
             py::arg("position") = ThreeVector(), py::arg("status") = 0,
             py::arg("mothers") = std::vector<int>(),
             py::arg("daughters") = std::vector<int>())
        // Setters
        .def("set_position", &Particle::SetPosition)
        .def("set_momentum", &Particle::SetMomentum)
        .def("set_status", &Particle::SetStatus)
        .def("set_mothers", &Particle::SetMothers)
        .def("set_daughters", &Particle::SetDaughters)
        .def("add_mother", &Particle::AddMother)
        .def("add_daughter", &Particle::AddDaughter)
        .def("set_formation_zone", &Particle::SetFormationZone)
        // Getters
        .def("pid", &Particle::ID)
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
        .def("get_distance_traveled", &Particle::GetDistanceTraveled)
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
