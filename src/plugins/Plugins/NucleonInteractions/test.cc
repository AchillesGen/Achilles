#include "Achilles/Interactions.hh"
#include "Achilles/Particle.hh"
#include "Plugins/InteractionLoader.hh"
#include "spdlog/spdlog.h"

int main() {
    achilles::InteractionLoader::LoadInteractions({"."});
    spdlog::info("Finished Registering");

    std::shared_ptr<achilles::Interactions> tmp =
        achilles::InteractionFactory::Create("ConstantInteractions", "10");
    achilles::Particle p1{};

    spdlog::info("Constant Cross-Section = {}", tmp->CrossSection(p1, p1));

    return 0;
}
