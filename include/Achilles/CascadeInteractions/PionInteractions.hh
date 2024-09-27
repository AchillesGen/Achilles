#ifndef ACHILLES_CASCADEINTERACTIONS_PIONINTERACTIONS
#define ACHILLES_CASCADEINTERACTIONS_PIONINTERACTIONS

#include "Achilles/Interactions.hh"

namespace achilles {

class PionInteraction : public Interaction, RegistrableInteraction<PionInteraction> {
  public:
    ///@name Constructors and Destructors
    ///@{

    /// Initialize PionAbsorptionTwoStep class
    PionInteraction() = default;
    PionInteraction(const YAML::Node &);

    /// Generate a object. This is used in the InteractionFactory.
    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<PionInteraction>(data);
    }

    /// Default Destructor
    ~PionInteraction() override = default;
    ///@}

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return PionInteraction::Name(); }
    static std::string Name() { return "PionInteraction"; }

    // List off all InitialStates
    std::vector<std::pair<PID, PID>> InitialStates() const override;

    InteractionResults CrossSection(Event &, size_t, size_t) const override;

    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

  private:
    std::unique_ptr<Interaction> hard_scatter;
    // TODO: Determine if this can remain an Interaction ptr or needs to be an Absorption ptr
    std::unique_ptr<Interaction> absorption;
};

} // namespace achilles

#endif
