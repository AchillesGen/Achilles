#include "Achilles/PyBindings.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ThreeVector.hh"

#include "plugins/InteractionLoader.hh"

namespace py = pybind11;
using namespace achilles;

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
    static std::unique_ptr<Interactions> Create(const YAML::Node&) {
        return std::make_unique<PyInteraction>(); 
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
    ThreeVector MakeMomentum(bool samePID,
            const double& pcm, const std::array<double, 2>& rans) const override {
        PYBIND11_OVERLOAD_PURE(
            ThreeVector,                // Return type
            Interactions,               // Parent class
            MakeMomentum,               // Name of function in C++ (must match Python name)
            samePID, pcm, rans    // Argument(s)
        );
    }

private:
    static bool registered;
};

REGISTER_INTERACTION(PyInteraction);

void InteractionsModule(py::module &m) {
    py::class_<Interactions, PyInteraction,
               std::shared_ptr<Interactions>>(m, "Interactions")
        .def(py::init<>())
        .def_static("create", [](const YAML::Node& data){
                return InteractionFactory::Create(data);})
        .def("cross_section", &Interactions::CrossSection)
        .def("make_momentum", &Interactions::MakeMomentum);

    py::class_<InteractionFactory>(m, "InteractionFactory")
        .def_static("register", &InteractionFactory::Register)
        .def_static("create", &InteractionFactory::Create)
        .def_static("list", &InteractionFactory::ListInteractions);

    py::class_<InteractionLoader>(m, "InteractionLoader")
        .def_static("load_plugins", &InteractionLoader::LoadInteractions);

}
