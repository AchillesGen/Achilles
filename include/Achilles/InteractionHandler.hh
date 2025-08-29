#pragma once

#include <functional>
#include <map>
#include <set>
#include <utility>

#include "Achilles/Interactions.hh"
#include "Achilles/Utilities.hh"
#include "yaml-cpp/yaml.h"

namespace achilles {

class Particle;
class PID;

class InteractionHandler {
    using Results = std::vector<InteractionResult>;
    using Channel = std::pair<PID, PID>;

  public:
    InteractionHandler() = default;
    InteractionHandler(const InteractionHandler &) = delete;
    InteractionHandler &operator=(const InteractionHandler &) = delete;
    InteractionHandler(InteractionHandler &&) = default;
    InteractionHandler &operator=(InteractionHandler &&) = default;
    double TotalCrossSection(Event &, size_t, size_t) const;
    std::vector<InteractionResult> CrossSection(Event &, size_t, size_t) const;
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
            // TODO: Move this logic to be better contained in the Settings or Factory classes
            try {
                auto interaction_obj = achilles::InteractionFactory::Initialize(
                    interaction["Name"].as<std::string>(), interaction["Options"]);
                handler.RegisterInteraction(std::move(interaction_obj));
            } catch(std::out_of_range &e) {
                spdlog::error(
                    "InteractionHandler: Requested interaction \"{}\", did you mean \"{}\"",
                    interaction["Name"].as<std::string>(),
                    achilles::GetSuggestion(achilles::InteractionFactory::List(),
                                            interaction["Name"].as<std::string>()));
                spdlog::error(
                    "InteractionHandler: Run `achilles --display-int-models` to see all options");
                exit(-1);
            }
        }
        return true;
    }
};
} // namespace YAML
