#pragma once

#include <functional>
#include <map>
#include <set>
#include <utility>

#include "Achilles/Interactions.hh"
#include "yaml-cpp/yaml.h"

namespace achilles {

class Particle;
class PID;

struct pid_compare {
    bool operator()(const std::pair<PID, PID> &lhs, const std::pair<PID, PID> &rhs) const;
};

class InteractionHandler {
    using Results = std::vector<InteractionResult>;
    using Channel = std::pair<PID, PID>;

  public:
    InteractionHandler() = default;
    InteractionHandler(const InteractionHandler &) = delete;
    InteractionHandler &operator=(const InteractionHandler &) = delete;
    InteractionHandler(InteractionHandler &&) = default;
    InteractionHandler &operator=(InteractionHandler &&) = default;
    double TotalCrossSection(const Particle &, const Particle &) const;
    std::vector<InteractionResult> CrossSection(const Particle &, const Particle &) const;
    std::vector<Particle> GenerateMomentum(const Particle &, const Particle &,
                                           const std::vector<PID> &, Random &&) const;
    std::set<Channel, pid_compare> RegisteredInteractions() const {
        return m_registered_interactions;
    }
    bool EnabledChannel(const PID &p1, const PID &p2) const;
    void RegisterInteraction(std::unique_ptr<Interaction>);
    InteractionResult SelectChannel(const std::vector<InteractionResult> &, double) const;

  private:
    std::vector<std::unique_ptr<Interaction>> m_interactions;
    std::set<Channel, pid_compare> m_registered_interactions;
    std::map<Channel, size_t, pid_compare> m_interaction_indices;
};

} // namespace achilles

namespace YAML {
template <> struct convert<achilles::InteractionHandler> {
    static bool decode(const Node &node, achilles::InteractionHandler &handler) {
        handler = achilles::InteractionHandler();
        for(const auto &interaction : node) {
            auto interaction_obj = achilles::InteractionFactory::Initialize(
                interaction["Name"].as<std::string>(), interaction["Options"]);
            handler.RegisterInteraction(std::move(interaction_obj));
        }
        return true;
    }
};
} // namespace YAML
