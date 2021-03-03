#include "nuchic/InteractionComponent.hh"
#include "spdlog/spdlog.h"

using namespace nuchic;

bool InteractionComponentFactory::Register(const std::string& name,
                                  InteractionComponentFactory::TCreateMethod funcCreate) {
    auto it = methods().find(name);
    if(it == methods().end()) {
        methods()[name] = funcCreate;
        spdlog::info("Registered: {}", name);
        return true;
    }
    return false;
}

std::unique_ptr<InteractionComponent> InteractionComponentFactory::Create(const YAML::Node& node) {
    auto name = node["Name"].as<std::string>();
    auto it = methods().find(name);
    if(it != methods().end())
        return it -> second(node);

    return nullptr;
}

void InteractionComponentFactory::ListInteractions() {
    fmt::print("+{:-^30s}+\n|{:^30s}|\n+{:-^30s}+\n", "", "Interactions", "");
    for(auto method : methods()) {
        fmt::print(" - {:30s}\n", method.first);
    }
}
