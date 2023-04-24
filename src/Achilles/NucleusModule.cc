#include "Achilles/FourVector.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/PyBindings.hh"
#include "Achilles/ThreeVector.hh"

// TODO: Deal with creating Nucleus in python since pybind11 does not like
// unique_ptr as arguments
void NucleusModule(py::module &m) {
    py::class_<achilles::Nucleus, std::shared_ptr<achilles::Nucleus>> nucleus(m, "Nucleus",
                                                                              py::module_local());
    // Constructors
    nucleus
        // Setters
        .def("set_nucleons", &achilles::Nucleus::SetNucleons)
        .def("set_binding_energy", &achilles::Nucleus::SetBindingEnergy)
        .def("set_fermi_momentum", &achilles::Nucleus::SetFermiMomentum)
        .def("set_potential", &achilles::Nucleus::SetPotential)
        // .def("set_density", &achilles::Nucleus::SetDensity)
        .def("set_radius", &achilles::Nucleus::SetRadius)
        // Getters
        .def("nucleons", &achilles::Nucleus::Nucleons)
        .def("protons", &achilles::Nucleus::Protons)
        .def("neutrons", &achilles::Nucleus::Neutrons)
        .def("n_nucleons", &achilles::Nucleus::NNucleons)
        .def("n_protons", &achilles::Nucleus::NProtons)
        .def("n_neutrons", &achilles::Nucleus::NNeutrons)
        .def("binding_energy", &achilles::Nucleus::BindingEnergy)
        .def("fermi_momentum", &achilles::Nucleus::FermiMomentum, py::arg("position") = 0.0)
        .def("potential_energy", &achilles::Nucleus::PotentialEnergy)
        .def("radius", &achilles::Nucleus::Radius)
        // Functions
        .def("escape", &achilles::Nucleus::Escape)
        .def("generate_config", &achilles::Nucleus::GenerateConfig)
        .def("generate_momentum", &achilles::Nucleus::GenerateMomentum, py::arg("position") = 0.0)
        // String Methods
        .def("__str__", &achilles::Nucleus::ToString)
        .def("__repr__", &achilles::Nucleus::ToString);
    // Static Methods
    // .def_static("make_nucleus", &achilles::Nucleus::MakeNucleus);

    py::enum_<achilles::Nucleus::FermiGasType>(nucleus, "FermiGas", py::module_local())
        .value("Local", achilles::Nucleus::FermiGasType::Local)
        .value("Global", achilles::Nucleus::FermiGasType::Global)
        .export_values();
}
