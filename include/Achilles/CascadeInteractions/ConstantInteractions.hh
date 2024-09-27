#ifndef ACHILLES_CASCADEINTERACTIONS_CONSTANTINTERACTIONS
#define ACHILLES_CASCADEINTERACTIONS_CONSTANTINTERACTIONS

#include <map>

#include "Achilles/Interactions.hh"
#include "Achilles/ParticleInfo.hh"

namespace achilles {

class ThreeVector;

/// Class for implementing an interaction model based on the Geant4 cross-section data. This
/// interaction model contains information about the angular distribution of pp, pn, and nn
/// interactions that occur during the intranuclear cascade.
class ConstantInteraction : public Interaction, RegistrableInteraction<ConstantInteraction> {
  public:
    ///@name Constructors and Destructors
    ///@{

    /// Initialize ConstantInteractions class for testing purposes.
    ConstantInteraction() = default;
    /// Initialize ConstantInteractions class. This loads data from an input file
    ///@param node: The constant interactions to load
    ConstantInteraction(const YAML::Node &);
    ConstantInteraction(const ConstantInteraction &) = default;
    ConstantInteraction(ConstantInteraction &&) = default;
    ConstantInteraction &operator=(const ConstantInteraction &) = default;
    ConstantInteraction &operator=(ConstantInteraction &&) = default;

    /// Generate a ConstantInteractions object. This is used in the InteractionFactory.
    ///@param data: The location of the data file to load containing the Geant4 cross-sections
    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<ConstantInteraction>(data);
    }

    /// Default Destructor
    ~ConstantInteraction() override = default;
    ///@}

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return ConstantInteraction::Name(); }
    static std::string Name() { return "ConstantInteraction"; }

    // These functions are defined in the base class
    std::vector<std::pair<PID, PID>> InitialStates() const override { return m_states; }
    InteractionResults CrossSection(Event &, size_t, size_t) const override;
    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

    // These functions are for testing only
    void AddInteraction(const std::pair<PID, PID> &state, const InteractionResults &results) {
        if(m_interactions.find(state) != m_interactions.end()) {
            auto msg =
                fmt::format("Initial state: [{}, {}] already exists", state.first, state.second);
            throw std::runtime_error("Initial state already exists");
        }
        m_interactions[state] = results;
        m_states.push_back(state);
    }

  private:
    ThreeVector MakeMomentum(bool, double, const std::vector<double> &) const;

    // Variables
    std::map<std::pair<PID, PID>, InteractionResults, pid_compare> m_interactions;
    std::vector<std::pair<PID, PID>> m_states;
};

} // namespace achilles

#endif
