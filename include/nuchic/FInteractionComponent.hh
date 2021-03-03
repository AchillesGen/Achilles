#ifndef FINTERACTION_COMPONENT_HH
#define FINTERACTION_COMPONENT_HH

#include "nuchic/InteractionComponent.hh"
#include "nuchic/Particle.hh"

extern "C" {
    using namespace nuchic;

    // Procedures implemented in Fortran to interface with Interactions
    void InitializeInteraction(const char*);
    double CrossSectionFortran(const Particle*, const Particle*);
    void CrossSectionsFortran(const Particle*, const Particle*, double**, int*);
    void GenerateFinalStateFortran(const Particle*, const Particle*, Particle**, int*);
}

namespace nuchic {

class FortranInteractionComponent : public InteractionComponent {
    public:
        FortranInteractionComponent(const YAML::Node &node) {
            auto name = node["Name"].as<std::string>();
            size_t len = name.size();
            auto cname = std::unique_ptr<char>(new char[len]);
            strcpy(cname.get(), name.c_str());
            InitializeInteraction(cname.get());
        };

        static std::unique_ptr<InteractionComponent> Create(const YAML::Node &data) {
            return std::make_unique<FortranInteractionComponent>(data);
        }

        static std::string GetName() { return "FortranInteraction"; }
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle &part1, const Particle &part2) const override {
            return CrossSectionFortran(&part1, &part2);
        }
        std::vector<double> CrossSections(const Particle &part1, const Particle &part2) const override {
            double *results = nullptr;
            int size{};
            CrossSectionsFortran(&part1, &part2, &results, &size);
            return std::vector<double>(results, results+size);
        }
        std::vector<Particle> GenerateFinalState(const Particle &part1, const Particle &part2) const override {
            Particle *results = nullptr;
            int size{};
            GenerateFinalStateFortran(&part1, &part2, &results, &size); 
            return std::vector<Particle>(results, results+size);
        }

        InteractionComponentType InteractionType() const override { return InteractionComponentType::NucleonNucleon; }

    private:
        static bool registered;
};

}

#endif
