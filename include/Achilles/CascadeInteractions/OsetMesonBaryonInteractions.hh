#ifndef ACHILLES_CASCADEINTERACTIONS_OsetMESONBARYONINTERACTIONS
#define ACHILLES_CASCADEINTERACTIONS_OsetMESONBARYONINTERACTIONS

#include "Achilles/Interactions.hh"
#include "Achilles/OsetCrossSections.hh"
#include "Achilles/Particle.hh"

namespace achilles {

class OsetMesonBaryonInteraction : public Interaction, RegistrableInteraction<OsetMesonBaryonInteraction> {
  public:
    ///@name Constructors and Destructors
    ///@{

    /// Initialize OsetMesonBaryonInteraction class.
    OsetMesonBaryonInteraction() = default;
    OsetMesonBaryonInteraction(const YAML::Node &);

    /// Generate a object. This is used in the InteractionFactory.
    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<OsetMesonBaryonInteraction>(data);
    }

    /// Default Destructor
    ~OsetMesonBaryonInteraction() override = default;
    ///@}

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return OsetMesonBaryonInteraction::Name(); }
    static std::string Name() { return "OsetMesonBaryonInteraction"; }

    // List off all InitialStates
    std::vector<std::pair<PID, PID>> InitialStates() const override;

    // Cross sections and pairs of PIDs for all possible final states
    InteractionResults CrossSection(Event &, size_t, size_t) const override;

    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

  private:
    // Variables
    std::map<std::pair<PID,PID>,std::vector<std::pair<PID,PID>>> in_out_states;

  protected:
    OsetCrossSection Oset; // Contains Oset absorption cross sections
};

} // namespace achilles

#endif
