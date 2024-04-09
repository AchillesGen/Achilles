#pragma once

#include <functional>
#include <map>
#include <set>
#include <utility>

namespace achilles {

class Particle;
class PID;
class Interaction;
struct InteractionResult;

struct pid_compare {
    bool operator()(const std::pair<PID, PID> &lhs, const std::pair<PID, PID> &rhs) const;
};

class InteractionHandler {
    using Results = std::vector<InteractionResult>;
    using Channel = std::pair<PID, PID>;
    using CrossSectionFunction = std::function<Results(const Particle &, const Particle &)>;
    using PhaseSpaceFunction = std::function<std::vector<Particle>(
        const Particle &, const Particle &, const std::vector<PID> &, const std::vector<double> &)>;

  public:
    std::vector<InteractionResult> CrossSection(const Particle &, const Particle &) const;
    std::vector<Particle> GenerateMomentum(const Particle &, const Particle &,
                                           const std::vector<PID> &,
                                           const std::vector<double> &) const;
    std::set<Channel, pid_compare> RegisteredInteractions() const {
        return m_registered_interactions;
    }
    bool EnabledChannel(const PID &p1, const PID &p2) const;
    void RegisterInteraction(const Interaction &);
    InteractionResult SelectChannel(const std::vector<InteractionResult> &, double) const;

  private:
    void RegisterCrossSection(const PID &, const PID &, CrossSectionFunction);
    void RegisterPhaseSpace(const PID &, const PID &, PhaseSpaceFunction);

    std::set<Channel, pid_compare> m_registered_interactions;
    std::map<Channel, CrossSectionFunction, pid_compare> m_cross_sections;
    std::map<Channel, PhaseSpaceFunction, pid_compare> m_phase_space;
};

} // namespace achilles
