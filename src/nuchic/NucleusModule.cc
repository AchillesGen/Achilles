#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/functional.h"

#include "nuchic/FourVector.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"

namespace py = pybind11;

PYBIND11_MODULE(nucleus, m) {
    py::object vectors = static_cast<py::object>(py::module::import("vectors"));
    py::object particle = static_cast<py::object>(py::module::import("particle"));

    py::class_<nuchic::Nucleus, std::shared_ptr<nuchic::Nucleus>>(m, "Nucleus")
        // Constructors
        .def(py::init<const int&, const int&, const double&,
                      const double&, const std::function<nuchic::Particles()>&>(),
                      py::arg("Z"), py::arg("A"), py::arg("binding"),
                      py::arg("kf"), py::arg("density") = std::function<nuchic::Particles()>())
        // Setters
        .def("set_nucleons", &nuchic::Nucleus::SetNucleons)
        .def("set_binding_energy", &nuchic::Nucleus::SetBindingEnergy)
        .def("set_fermi_momentum", &nuchic::Nucleus::SetFermiMomentum)
        .def("set_potential", &nuchic::Nucleus::SetPotential)
        .def("set_density", &nuchic::Nucleus::SetDensity)
        .def("set_radius", &nuchic::Nucleus::SetRadius)
        // Getters
        .def("nucleons", &nuchic::Nucleus::Nucleons)
        .def("protons", &nuchic::Nucleus::Protons)
        .def("neutrons", &nuchic::Nucleus::Neutrons)
        .def("n_nucleons", &nuchic::Nucleus::NNucleons)
        .def("n_protons", &nuchic::Nucleus::NProtons)
        .def("n_neutrons", &nuchic::Nucleus::NNeutrons)
        .def("binding_energy", &nuchic::Nucleus::BindingEnergy)
        .def("fermi_momentum", &nuchic::Nucleus::FermiMomentum)
        .def("potential_energy", &nuchic::Nucleus::PotentialEnergy)
        .def("radius", &nuchic::Nucleus::Radius)
        // Functions
        .def("escape", &nuchic::Nucleus::Escape)
        .def("generate_config", &nuchic::Nucleus::GenerateConfig)
        .def("generate_momentum", &nuchic::Nucleus::GenerateMomentum)
        // String Methods
        .def("__str__", &nuchic::Nucleus::ToString)
        .def("__repr__", &nuchic::Nucleus::ToString)
        // Static Methods
        .def_static("make_nucleus", &nuchic::Nucleus::MakeNucleus);
}

