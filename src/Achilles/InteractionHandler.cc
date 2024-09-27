#include "Achilles/InteractionHandler.hh"
#include "Achilles/Event.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Random.hh"

#include <fmt/format.h>

using achilles::InteractionHandler;
using achilles::InteractionResult;
using achilles::Particle;

double InteractionHandler::TotalCrossSection(Event &event, size_t part1, size_t part2) const {
    const auto &p1 = event.Hadrons()[part1];
    const auto &p2 = event.Hadrons()[part2];

    auto key = std::make_pair(p1.ID(), p2.ID());
    auto it = m_interaction_indices.find(key);
    if(it == m_interaction_indices.end()) {
        auto msg =
            fmt::format("Cross section not registered for particles: [{}, {}]", p1.ID(), p2.ID());
        throw std::runtime_error(msg);
    }
    return m_interactions[it->second]->TotalCrossSection(event, part1, part2);
}

std::vector<InteractionResult> InteractionHandler::CrossSection(Event &event, size_t part1,
                                                                size_t part2) const {
    const auto &p1 = event.Hadrons()[part1];
    const auto &p2 = event.Hadrons()[part2];

    auto key = std::make_pair(p1.ID(), p2.ID());
    auto it = m_interaction_indices.find(key);
    if(it == m_interaction_indices.end()) {
        auto msg =
            fmt::format("Cross section not registered for particles: [{}, {}]", p1.ID(), p2.ID());
        throw std::runtime_error(msg);
    }
    return m_interactions[it->second]->CrossSection(event, part1, part2);
}

std::vector<Particle> InteractionHandler::GenerateMomentum(const Particle &p1, const Particle &p2,
                                                           const std::vector<PID> &out_pids,
                                                           Random &&random) const {
    auto key = std::make_pair(p1.ID(), p2.ID());
    auto it = m_interaction_indices.find(key);
    if(it == m_interaction_indices.end()) {
        auto msg =
            fmt::format("Phase space not registered for particles: [{}, {}]", p1.ID(), p2.ID());
        throw std::runtime_error(msg);
    }
    return m_interactions[it->second]->GenerateMomentum(p1, p2, out_pids, random);
}

bool InteractionHandler::EnabledChannel(const PID &p1, const PID &p2) const {
    return m_registered_interactions.find(std::make_pair(p1, p2)) !=
           m_registered_interactions.end();
}

void InteractionHandler::RegisterInteraction(std::unique_ptr<Interaction> interaction) {
    size_t idx = m_interactions.size();
    m_interactions.push_back(std::move(interaction));
    for(const auto &pid : m_interactions.back()->InitialStates()) {
        if(EnabledChannel(pid.first, pid.second)) {
            throw std::runtime_error(fmt::format(
                "Interaction already registered for particles: [{}, {}]", pid.first, pid.second));
        }

        spdlog::debug("InteractionHandler: Registering interaction for particles: [{}, {}]",
                      pid.first, pid.second);
        m_registered_interactions.insert(pid);
        m_interaction_indices[pid] = idx;
    }
}

InteractionResult InteractionHandler::SelectChannel(const std::vector<InteractionResult> &channels,
                                                    double rand) const {
    if(channels.empty()) { throw std::runtime_error("No channels to select from"); }

    double total = 0.0;
    for(const auto &channel : channels) { total += channel.cross_section; }

    for(const auto &channel : channels) {
        rand -= channel.cross_section / total;
        if(rand <= 0.0) { return channel; }
    }

    throw std::runtime_error("No channel selected");
}
