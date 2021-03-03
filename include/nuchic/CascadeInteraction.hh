#ifndef CASCADE_INTERACTION_HH
#define CASCADE_INTERACTION_HH

#include "nuchic/InteractionComponentEnums.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace nuchic {

class ParticleInfo;
class InteractionComponent;

class CascadeInteraction {
    public:
        CascadeInteraction() = default;
        CascadeInteraction(const CascadeInteraction&) = delete;
        CascadeInteraction(CascadeInteraction&&);
        CascadeInteraction& operator=(const CascadeInteraction&) = delete;
        CascadeInteraction& operator=(CascadeInteraction&&) = default;
        ~CascadeInteraction() = default;

        void AddComponent(const YAML::Node&);
        InteractionComponent* GetComponent(nuchic::ParticleInfo) const;

        // For testing
        void AddComponent(std::unique_ptr<InteractionComponent>);
    private:
        std::map<InteractionComponentType, std::unique_ptr<InteractionComponent>> components;
};

}

namespace YAML {
template<>
struct convert<nuchic::CascadeInteraction> {
    static bool decode(const Node &components, nuchic::CascadeInteraction &cinteraction) {
        for(const auto& component : components) {
            cinteraction.AddComponent(component);
        }

        return true;
    }
};
}

#endif
