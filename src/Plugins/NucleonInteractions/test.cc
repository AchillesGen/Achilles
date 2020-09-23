#include "nuchic/Interactions.hh"
#include "Plugins/InteractionLoader.hh"
#include "spdlog/spdlog.h"
#include "nuchic/Particle.hh"

int main() {
    nuchic::InteractionLoader::LoadInteractions({"."}); 
    spdlog::info("Finished Registering");

    std::shared_ptr<nuchic::Interactions> tmp = nuchic::InteractionFactory::Create("ConstantInteractions", "10");
    nuchic::Particle p1{};
    
    spdlog::info("Constant Cross-Section = {}", tmp->CrossSection(p1, p1));

    return 0;
}
