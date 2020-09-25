#include "nuchic/Interactions.hh"
#include "spdlog/spdlog.h"

using namespace nuchic;

bool InteractionFactory::Register(const std::string& name,
                                  InteractionFactory::TCreateMethod funcCreate) {
    auto it = methods().find(name);
    if(it == methods().end()) {
        methods()[name] = funcCreate;
        if(name != "PyWrapper")
            spdlog::info("Registered: {}", name);
        return true;
    }
    return false;
}

std::shared_ptr<Interactions> InteractionFactory::Create(const std::string& name,
                                                         const std::string& filename) {
    if(name.substr(0, 2) == "Py") {
        auto it = methods().find("PyWrapper");
        if(it != methods().end())
            return it -> second(name.substr(2));
    }
    auto it = methods().find(name);
    if(it != methods().end())
        return it -> second(filename);

    return nullptr;
}

void InteractionFactory::ListInteractions() {
    fmt::print("+{:-^30s}+\n|{:^30s}|\n+{:-^30s}+\n", "", "Interactions", "");
    fmt::print("C++ Interactions:\n");
    fmt::print("-----------------\n");
    for(auto method : methods()) {
        if(method.first != "PyWrapper")
            fmt::print(" - {:30s}\n", method.first);
    }
    fmt::print("\n");
}
