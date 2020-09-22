#include "nuchic/PyBindings.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"

void NucleusModule(py::module &m) {
    py::class_<nuchic::Nucleus, std::shared_ptr<nuchic::Nucleus>> nucleus(m, "Nucleus", py::module_local());
    // Constructors
    nucleus.def(py::init<const std::size_t&, const std::size_t&, const double&, const double&,
                const std::string&, const nuchic::Nucleus::FermiGasType&,
                const std::function<nuchic::Particles()>&>(),
                py::arg("Z"), py::arg("A"), py::arg("binding"), py::arg("kf"),
                py::arg("density_file"),
                py::arg("fg_type"),
                py::arg("density") = std::function<nuchic::Particles()>())
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
        .def("fermi_momentum", &nuchic::Nucleus::FermiMomentum,
             py::arg("position") = 0.0)
        .def("potential_energy", &nuchic::Nucleus::PotentialEnergy)
        .def("radius", &nuchic::Nucleus::Radius)
        // Functions
        .def("escape", &nuchic::Nucleus::Escape)
        .def("generate_config", &nuchic::Nucleus::GenerateConfig)
        .def("generate_momentum", &nuchic::Nucleus::GenerateMomentum,
             py::arg("position") = 0.0)
        // String Methods
        .def("__str__", &nuchic::Nucleus::ToString)
        .def("__repr__", &nuchic::Nucleus::ToString)
        // Static Methods
        .def_static("make_nucleus", &nuchic::Nucleus::MakeNucleus);

    py::enum_<nuchic::Nucleus::FermiGasType>(nucleus, "FermiGas", py::module_local())
        .value("Local", nuchic::Nucleus::FermiGasType::Local)
        .value("Global", nuchic::Nucleus::FermiGasType::Global)
        .export_values();
}

