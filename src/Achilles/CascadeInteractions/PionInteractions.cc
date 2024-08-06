#include "Achilles/CascadeInteractions/PionInteractions.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"

using namespace achilles;

PionInteraction::PionInteraction(const YAML::Node &node) {
    auto hard_scatter_node = node["HardScatter"];
    auto absorption_node = node["Absorption"];
    hard_scatter = InteractionFactory::Initialize(hard_scatter_node["Name"].as<std::string>(),
                                                  hard_scatter_node["Options"]);
    absorption = InteractionFactory::Initialize(absorption_node["Name"].as<std::string>(),
                                                absorption_node["Options"]);
}

std::vector<std::pair<PID, PID>> PionInteraction::InitialStates() const {
    return {{PID::pion0(), PID::proton()},  {PID::pionp(), PID::proton()},
            {-PID::pionp(), PID::proton()}, {PID::pion0(), PID::neutron()},
            {PID::pionp(), PID::neutron()}, {-PID::pionp(), PID::neutron()}};
}

std::vector<Particle> PionInteraction::GenerateMomentum(const Particle &particle1,
                                                        const Particle &particle2,
                                                        const std::vector<PID> &out_pids,
                                                        Random &random) const {
    bool is_meson = false;
    for(const auto &pid : out_pids) {
        if(ParticleInfo(pid).IsBoson()) {
            is_meson = true;
            break;
        }
    }

    if(is_meson) return hard_scatter->GenerateMomentum(particle1, particle2, out_pids, random);
    return absorption->GenerateMomentum(particle1, particle2, out_pids, random);
}

InteractionResults PionInteraction::CrossSection(Event &event, size_t part1, size_t part2) const {
    auto Quasielastic_xsec = hard_scatter->CrossSection(event, part1, part2);
    auto Absorption_xsec = absorption->CrossSection(event, part1, part2);

    InteractionResults results = Quasielastic_xsec;
    results.insert(results.begin(), Absorption_xsec.begin(), Absorption_xsec.end());

    return results;
}
