#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include "pybind11/functional.h"

#include "nuchic/Interactions.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"

#include "Plugins/InteractionLoader.hh"

namespace py = pybind11;
using namespace nuchic;

class PyWrapper : public Interactions {
public:
    // Inherit the constructors 
    // using Interactions::Interactions;
    PyWrapper() = default;
    PyWrapper(const PyWrapper&) = default;
    PyWrapper(PyWrapper&&) = default;
    PyWrapper(py::object pyClass) : m_pyClass(std::move(pyClass)) {}
    PyWrapper& operator=(const PyWrapper&) = default;
    PyWrapper& operator=(PyWrapper&&) = default;
    /// Default Destructor
    ~PyWrapper() override = default;

    static bool IsRegistered() noexcept { return registered; }
    static std::string GetName() { return "PyWrapper"; }
    static std::unique_ptr<Interactions> Create(const std::string &name) {
        py::module module = py::module::import("nuchic.pyinteractions");
        py::object obj = module.attr("PyInteraction")(name, 10);
        return std::make_unique<PyWrapper>(obj); 
    }

    // Trampoline for CrossSection
    double CrossSection(const Particle& part1, const Particle& part2) const override {
        return m_pyClass.attr("CrossSection")(part1, part2).cast<double>();
    }

    // Trampoline for MakeMomentum
    ThreeVector MakeMomentum(bool samePID, const double& p1CM,
            const double& pcm, const std::array<double, 2>& rans) const override {
        return m_pyClass.attr("MakeMomentum")(samePID, p1CM, pcm, rans).cast<ThreeVector>();
    }

private:
    py::object m_pyClass;
    static bool registered;
};


REGISTER_INTERACTION(PyWrapper);

class PyInteraction : public Interactions {
public:
    // Inherit the constructors 
    // using Interactions::Interactions;
    PyInteraction() = default;
    PyInteraction(const PyInteraction&) = default;
    PyInteraction(PyInteraction&&) = default;
    PyInteraction& operator=(const PyInteraction&) = default;
    PyInteraction& operator=(PyInteraction&&) = default;
    /// Default Destructor
    ~PyInteraction() override = default;

    static bool IsRegistered() noexcept { return registered; }
    static std::string GetName() { return "PyInteraction"; }
    static std::unique_ptr<Interactions> Create(const std::string &name) {
        py::module module = py::module::import("nuchic.run_modes");
        std::cout << "Name: " << name << std::endl;
        py::object obj = module.attr("PyInteraction")(name, 10);
        std::cout << "Making unique" << std::endl;
        return std::make_unique<PyInteraction>(obj.cast<PyInteraction>()); 
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

private:
    static bool registered;
};

PYBIND11_MODULE(interactions, m) {
    py::class_<Interactions, PyWrapper,
               std::shared_ptr<Interactions>>(m, "Interactions")
        .def(py::init<>())
        .def_static("create", [](const std::string& name, const std::string& data){
                return InteractionFactory::Create(name, data);})
        .def("cross_section", &Interactions::CrossSection)
        .def("make_momentum", &Interactions::MakeMomentum);

    py::class_<InteractionFactory>(m, "InteractionFactory")
        .def_static("register", &InteractionFactory::Register)
        .def_static("create", &InteractionFactory::Create)
        .def_static("list", &InteractionFactory::ListInteractions);

    py::class_<InteractionLoader>(m, "InteractionLoader")
        .def_static("load_plugins", &InteractionLoader::LoadInteractions);

}
