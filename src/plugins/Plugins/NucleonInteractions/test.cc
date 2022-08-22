#include "Achilles/Interactions.hh"
#include "Plugins/InteractionLoader.hh"
#include "spdlog/spdlog.h"
#include "Achilles/Particle.hh"

int main() {
    achilles::InteractionLoader::LoadInteractions({"."}); 
    spdlog::info("Finished Registering");

    std::shared_ptr<achilles::Interactions> tmp = achilles::InteractionFactory::Create("ConstantInteractions", "10");
    achilles::Particle p1{};
    
    spdlog::info("Constant Cross-Section = {}", tmp->CrossSection(p1, p1));

    return 0;
}
