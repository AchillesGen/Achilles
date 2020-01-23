#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include "pybind11/functional.h"

#include "nuchic/Interactions.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"

namespace py = pybind11;

class PyInteraction : public Interactions {
public:
    // Inherit the constructors 
    using Interactions::Interactions;

    // Trampoline for IsRegistered
    bool IsRegistered() const noexcept override {
        PYBIND11_OVERLOAD_PURE(
            bool,            // Return type
            Interactions,      // Parent class
            IsRegistered      // Name of function in C++ (must match Python name)
        );
    }

    // Trampoline for CrossSection
    double CrossSection(const Particle& part1, const Particle& part2) const override {
        PYBIND11_OVERLOAD_PURE(
            double,            // Return type
            Interactions,      // Parent class
            CrossSection,      // Name of function in C++ (must match Python name)
            part1, part2       // Argument(s)
        );
    }

    // Trampoline for MakeMomentum
    ThreeVector MakeMomentum(bool samePID, const double& p1CM,
            const double& pcm, const std::array<double, 2>& rans) const override {
        PYBIND11_OVERLOAD_PURE(
            ThreeVector,                // Return type
            Interactions,               // Parent class
            MakeMomentum,               // Name of function in C++ (must match Python name)
            samePID, p1CM, pcm, rans    // Argument(s)
        );
    }
};

PYBIND11_MODULE(interactions, m) {
    m.def("cross_section", &CrossSection);
    m.def("cross_section_lab", &CrossSectionLab);
    m.def("cross_section_angle", &CrossSectionAngle);
    m.def("make_momentum_angular", &MakeMomentumAngular);

    py::class_<Interactions, PyInteraction/*Trampoline*/, std::shared_ptr<Interactions>>(m, "Interactions")
        .def(py::init<>())
        .def_static("create", [](const std::string& name){return InteractionFactory::Instance().Create(name);})
        .def("CrossSection", &Interactions::CrossSection)
        .def("MakeMomentum", &Interactions::MakeMomentum);

    py::class_<InteractionFactory>(m, "InteractionFactory")
        .def_static("instance", &InteractionFactory::Instance)
        .def("register", &InteractionFactory::Register)
        .def("create", &InteractionFactory::Create);

    py::class_<GeantInteractions, Interactions, std::shared_ptr<GeantInteractions>>(m, "GeantInteractions")
        .def(py::init<const std::string&>())
        .def("CrossSection", &GeantInteractions::CrossSection)
        .def("MakeMomentum", &GeantInteractions::MakeMomentum);
}
