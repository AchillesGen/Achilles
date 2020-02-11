#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include "pybind11/functional.h"

#include "nuchic/Interactions.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"

namespace py = pybind11;

class PyInteraction : public nuchic::Interactions {
public:
    // Inherit the constructors 
    using nuchic::Interactions::Interactions;

    // Trampoline for IsRegistered
    bool IsRegistered() const noexcept override {
        PYBIND11_OVERLOAD_PURE(
            bool,            // Return type
            nuchic::Interactions,      // Parent class
            IsRegistered      // Name of function in C++ (must match Python name)
        );
    }

    // Trampoline for CrossSection
    double CrossSection(const nuchic::Particle& part1, const nuchic::Particle& part2) const override {
        PYBIND11_OVERLOAD_PURE(
            double,            // Return type
            nuchic::Interactions,      // Parent class
            CrossSection,      // Name of function in C++ (must match Python name)
            part1, part2       // Argument(s)
        );
    }

    // Trampoline for MakeMomentum
    nuchic::ThreeVector MakeMomentum(bool samePID, const double& p1CM,
            const double& pcm, const std::array<double, 2>& rans) const override {
        PYBIND11_OVERLOAD_PURE(
            nuchic::ThreeVector,                // Return type
            nuchic::Interactions,               // Parent class
            MakeMomentum,               // Name of function in C++ (must match Python name)
            samePID, p1CM, pcm, rans    // Argument(s)
        );
    }
};

PYBIND11_MODULE(interactions, m) {
    m.def("cross_section", &nuchic::CrossSection);
    m.def("cross_section_lab", &nuchic::CrossSectionLab);
    m.def("cross_section_angle", &nuchic::CrossSectionAngle);
    m.def("make_momentum_angular", &nuchic::MakeMomentumAngular);

    py::class_<nuchic::Interactions, PyInteraction,
               std::shared_ptr<nuchic::Interactions>>(m, "Interactions")
        .def(py::init<>())
        .def_static("create", [](const std::string& name, const std::string& data){
                return nuchic::InteractionFactory::Instance().Create(name, data);})
        .def("CrossSection", &nuchic::Interactions::CrossSection)
        .def("MakeMomentum", &nuchic::Interactions::MakeMomentum);

    py::class_<nuchic::InteractionFactory>(m, "InteractionFactory")
        .def_static("instance", &nuchic::InteractionFactory::Instance)
        .def("register", &nuchic::InteractionFactory::Register)
        .def("create", &nuchic::InteractionFactory::Create);
}
