#ifndef ACHILLES_CASCADEINTERACTIONS_MESONBARYONINTERACTIONS
#define ACHILLES_CASCADEINTERACTIONS_MESONBARYONINTERACTIONS

#include "Achilles/Interactions.hh"
#include "Achilles/MesonBaryonAmplitudes.hh"

namespace achilles {

class MesonBaryonInteraction : public Interaction, RegistrableInteraction<MesonBaryonInteraction> {
  public:
    ///@name Constructors and Destructors
    ///@{

    /// Initialize MesonBaryonInteraction class. This loads data from an input file
    MesonBaryonInteraction() = default;
    MesonBaryonInteraction(const YAML::Node &){};

    /// Generate a object. This is used in the InteractionFactory.
    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<MesonBaryonInteraction>(data);
    }

    /// Default Destructor
    ~MesonBaryonInteraction() override = default;
    ///@}

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return MesonBaryonInteraction::Name(); }
    static std::string Name() { return "MesonBaryonInteraction"; }

    // List off all InitialStates
    std::vector<std::pair<PID, PID>> InitialStates() const override;

    // Cross sections and pairs of PIDs for all possible final states
    InteractionResults CrossSection(Event &, size_t, size_t) const override;

    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

  private:
    // Variables
    MBAmplitudes m_amps; // Contains all the mesonBaryon information
};

} // namespace achilles

#endif
