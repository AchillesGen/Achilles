#include "nuchic/PyBindings.hh"
#include "nuchic/Interactions.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"

using namespace nuchic;

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
    static std::unique_ptr<Interactions> Create(const std::string&) {
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
    m.def("cross_section", &CrossSection);
    m.def("cross_section_lab", &CrossSectionLab);
    m.def("cross_section_angle", &CrossSectionAngle);
    m.def("make_momentum_angular", &MakeMomentumAngular);

    py::class_<Interactions, PyInteraction,
               std::shared_ptr<Interactions>>(m, "Interactions")
        .def(py::init<>())
        .def_static("create", [](const std::string& name, const std::string& data){
                return InteractionFactory::Create(name, data);})
        .def("cross_section", &Interactions::CrossSection)
        .def("make_momentum", &Interactions::MakeMomentum);

    py::class_<InteractionFactory>(m, "InteractionFactory")
        .def("register", &InteractionFactory::Register)
        .def("create", &InteractionFactory::Create)
        .def("list", &InteractionFactory::ListInteractions);
}
