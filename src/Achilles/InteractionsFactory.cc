#include "Achilles/Interactions.hh"
#include "spdlog/spdlog.h"

using namespace achilles;

bool InteractionFactory::Register(const std::string& name,
                                  InteractionFactory::TCreateMethod funcCreate) {
    auto it = methods().find(name);
    if(it == methods().end()) {
        methods()[name] = funcCreate;
        if(name != "PyInteraction")
            spdlog::info("Registered: {}", name);
        return true;
    }
    return false;
}

std::unique_ptr<Interactions> InteractionFactory::Create(const YAML::Node& node) {
    auto name = node["Name"].as<std::string>();
    auto it = methods().find(name);
    if(it != methods().end())
        return it -> second(node);

    spdlog::error("Interaction {} is undefined", name);
    throw std::runtime_error(fmt::format("Invalid Interaction Mode", name));
}

void InteractionFactory::Display() {
    fmt::print("+{:-^30s}+\n|{:^30s}|\n+{:-^30s}+\n", "", "Interactions", "");
    for(auto method : methods()) {
        fmt::print(" - {:30s}\n", method.first);
    }
}
