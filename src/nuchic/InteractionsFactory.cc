#include "nuchic/Interactions.hh"
#include "spdlog/spdlog.h"

using namespace nuchic;

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

std::shared_ptr<Interactions> InteractionFactory::Create(const std::string& name,
                                                         const std::string& filename) {
    auto it = methods().find(name);
    if(it != methods().end())
        return it -> second(filename);

    spdlog::error("Interaction {} is undefined", name);
    throw std::runtime_error(fmt::format("Invalid Interaction Mode", name));
}

void InteractionFactory::ListInteractions() {
    fmt::print("+{:-^30s}+\n|{:^30s}|\n+{:-^30s}+\n", "", "Interactions", "");
    for(auto method : methods()) {
        fmt::print(" - {:30s}\n", method.first);
    }
}
