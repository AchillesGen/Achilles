#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/functional.h"

#include "nuchic/FourVector.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"

namespace py = pybind11;

PYBIND11_MODULE(nucleus, m) {
    py::class_<Nucleus, std::shared_ptr<Nucleus>>(m, "Nucleus")
        // Constructors
        .def(py::init<const int&, const int&, const double&,
                      const double&, const std::function<Particles()>&>())
        // Setters
        .def("set_nucleons", &Nucleus::SetNucleons)
        .def("set_binding_energy", &Nucleus::SetBindingEnergy)
        .def("set_fermi_momentum", &Nucleus::SetFermiMomentum)
        .def("set_potential", &Nucleus::SetPotential)
        .def("set_density", &Nucleus::SetDensity)
        // Getters
        .def("nucleons", &Nucleus::Nucleons)
        .def("protons", &Nucleus::Protons)
        .def("neutrons", &Nucleus::Neutrons)
        .def("n_nucleons", &Nucleus::NNucleons)
        .def("n_protons", &Nucleus::NProtons)
        .def("n_neutrons", &Nucleus::NNeutrons)
        .def("binding_energy", &Nucleus::BindingEnergy)
        .def("fermi_momentum", &Nucleus::FermiMomentum)
        .def("potential_energy", &Nucleus::PotentialEnergy)
        // Functions
        .def("escape", &Nucleus::Escape)
        .def("generate_config", &Nucleus::GenerateConfig)
        .def("generate_momentum", &Nucleus::GenerateMomentum)
        // String Methods
        .def("__str__", &Nucleus::ToString)
        .def("__repr__", &Nucleus::ToString)
        // Static Methods
        .def_static("make_nucleus", &Nucleus::MakeNucleus);
}

