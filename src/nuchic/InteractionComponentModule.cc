#include "nuchic/PyBindings.hh"
#include "nuchic/InteractionComponent.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"

using namespace nuchic;

class PyInteraction : public InteractionComponent {
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
    static std::unique_ptr<InteractionComponent> Create(const YAML::Node&) {
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

    // Trampoline for CrossSections
    std::vector<double> CrossSections(const Particle &part1, const Particle &part2) const override {
        PYBIND11_OVERLOAD_PURE(
            std::vector<double>,            // Return type
            Interactions,                   // Parent class
            CrossSections,                  // Name of function in C++ (must match Python name)
            part1, part2                    // Argument(s)
        );
    }

    // Trampoline for CrossSections
    std::vector<Particle> GenerateFinalState(randutils::mt19937_rng &rng,
                                             const Particle &part1,
                                             const Particle &part2) const override {
        PYBIND11_OVERLOAD_PURE(
            std::vector<Particle>,            // Return type
            Interactions,                     // Parent class
            GenerateFinalState,               // Name of function in C++ (must match Python name)
            rng, part1, part2                 // Argument(s)
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
    m.def("cross_section", &CrossSection);
    m.def("cross_section_lab", &CrossSectionLab);
    m.def("cross_section_angle", &CrossSectionAngle);
    m.def("make_momentum_angular", &MakeMomentumAngular);

    py::class_<InteractionComponent, PyInteraction,
               std::shared_ptr<InteractionComponent>>(m, "Interactions")
        .def(py::init<>())
        .def_static("create", [](const YAML::Node& data){
                return InteractionComponentFactory::Create(data);})
        .def("cross_section", &InteractionComponent::CrossSection)
        .def("make_momentum", &InteractionComponent::MakeMomentum);

    py::class_<InteractionComponentFactory>(m, "InteractionFactory")
        .def("register", &InteractionComponentFactory::Register)
        .def("create", &InteractionComponentFactory::Create)
        .def("list", &InteractionComponentFactory::ListInteractions);
}
