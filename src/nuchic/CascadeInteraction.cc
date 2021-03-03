#include "nuchic/CascadeInteraction.hh"
#include "nuchic/InteractionComponent.hh"

nuchic::CascadeInteraction::CascadeInteraction(nuchic::CascadeInteraction &&other) {
    for(auto &pair : other.components) {
        components.insert({pair.first, std::move(pair.second)});
    }
}

void nuchic::CascadeInteraction::AddComponent(const YAML::Node &node) {
    auto component = nuchic::InteractionComponentFactory::Create(node);
    InteractionComponentType type = component->InteractionType();
    spdlog::debug("Adding component for {}", static_cast<int>(type));
    components.insert({type, std::move(component)});
}

void nuchic::CascadeInteraction::AddComponent(std::unique_ptr<InteractionComponent> component) {
    InteractionComponentType type = component -> InteractionType();
    components.insert({type, std::move(component)});
}

nuchic::InteractionComponent* nuchic::CascadeInteraction::GetComponent(nuchic::ParticleInfo info) const {
    nuchic::InteractionComponent *component = nullptr;
    if(info.IsNucleon()) {
        spdlog::debug("Using interaction component NucleonNucleon");
        component = components.at(InteractionComponentType::NucleonNucleon).get();
    } else if(info.IsPion()) {
        spdlog::debug("Using interaction component NucleonPion");
        component = components.at(InteractionComponentType::NucleonPion).get();
    }
    return component;
}
