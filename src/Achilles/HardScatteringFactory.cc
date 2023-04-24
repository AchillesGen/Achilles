#include "Achilles/HardScatteringFactory.hh"
#include "Achilles/HardScattering.hh"
#include "spdlog/spdlog.h"

using achilles::HardScattering;
using achilles::HardScatteringFactory;

bool HardScatteringFactory::Register(const std::string &name, TCreateMethod funcCreate) {
    auto it = methods().find(name);
    if(it == methods().end()) {
        methods()[name] = funcCreate;
        spdlog::info("Registered: {}", name);
        return true;
    }
    return false;
}

std::unique_ptr<HardScattering> HardScatteringFactory::Create(const std::string &name,
                                                              const YAML::Node &config,
                                                              achilles::RunMode mode) {
    auto it = methods().find(name);
    if(it != methods().end()) return it->second(config, mode);
    return nullptr;
}

void HardScatteringFactory::ListHardScattering() {
    fmt::print("+{:-^30s}+\n|{:^30s}|\n+{:-^30s}+\n", "", "HardScattering", "");
    for(auto method : methods()) { fmt::print(" - {:30s}\n", method.first); }
}
