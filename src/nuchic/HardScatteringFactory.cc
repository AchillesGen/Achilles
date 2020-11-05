#include "nuchic/HardScatteringFactory.hh"
#include "nuchic/HardScattering.hh"
#include "spdlog/spdlog.h"

using nuchic::HardScatteringFactory;
using nuchic::HardScattering;

bool HardScatteringFactory::Register(const std::string &name,
                                     TCreateMethod funcCreate) {
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
                                                              std::shared_ptr<Beam> beam,
                                                              std::shared_ptr<Nucleus> nucleus,
                                                              nuchic::RunMode mode) {
    auto it = methods().find(name);
    if(it != methods().end())
        return it -> second(config, beam, nucleus, mode);
    return nullptr;
}

void HardScatteringFactory::ListHardScattering() {
    fmt::print("+{:-^30s}+\n|{:^30s}|\n+{:-^30s}+\n", "", "HardScattering", "");
    for(auto method : methods()) {
        fmt::print(" - {:30s}\n", method.first);
    }
}
