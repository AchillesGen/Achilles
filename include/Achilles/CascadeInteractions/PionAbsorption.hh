#ifndef ACHILLES_CASCADEINTERACTIONS_PIONABSORPTION
#define ACHILLES_CASCADEINTERACTIONS_PIONABSORPTION

#include "Achilles/Interactions.hh"
#include "Achilles/OsetAbsorption.hh"
#include "Achilles/Particle.hh"

namespace achilles {

class PionAbsorption : public Interaction {
  public:
    /// Initialize PionAbsorption Helper class
    PionAbsorption() = default;
    PionAbsorption(const YAML::Node &);

    // List off all InitialStates
    std::vector<std::pair<PID, PID>> InitialStates() const override;

    InteractionResults CrossSection(Event &, size_t, size_t) const override;

  protected:
    virtual bool AllowedAbsorption(Event &, size_t, size_t) const = 0;
    mutable Particle absorption_partner;

    OsetAbsCrossSection Oset_abs; // Contains Oset absorption cross sections

    struct InteractionStates {
        PID absorption_partner;
        std::vector<PID> outgoing;
    };

    using InteractionMap = std::map<std::pair<PID, PID>, std::vector<InteractionStates>>;
    InteractionMap states;
};

class PionAbsorptionOneStep : public PionAbsorption, RegistrableInteraction<PionAbsorptionOneStep> {
  public:
    ///@name Constructors and Destructors
    ///@{

    /// Initialize PionAbsorptionOneStep class
    PionAbsorptionOneStep() = default;
    PionAbsorptionOneStep(const YAML::Node &){};

    /// Generate a object. This is used in the InteractionFactory.
    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<PionAbsorptionOneStep>(data);
    }

    /// Default Destructor
    ~PionAbsorptionOneStep() override = default;
    ///@}

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return PionAbsorptionOneStep::Name(); }
    static std::string Name() { return "PionAbsorptionOneStep"; }

    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

  private:
    std::pair<size_t, size_t> FindClosest(Event &, size_t, size_t) const;
    bool AllowedAbsorption(Event &, size_t, size_t) const override;
    mutable std::vector<std::pair<size_t, std::vector<PID>>> m_states;
};

class PionAbsorptionTwoStep : public PionAbsorption, RegistrableInteraction<PionAbsorptionTwoStep> {
  public:
    ///@name Constructors and Destructors
    ///@{

    /// Initialize PionAbsorptionTwoStep class
    PionAbsorptionTwoStep() = default;
    PionAbsorptionTwoStep(const YAML::Node &){};

    /// Generate a object. This is used in the InteractionFactory.
    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<PionAbsorptionTwoStep>(data);
    }

    /// Default Destructor
    ~PionAbsorptionTwoStep() override = default;
    ///@}

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return PionAbsorptionTwoStep::Name(); }
    static std::string Name() { return "PionAbsorptionTwoStep"; }

    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

  private:
    bool CheckHistory(Event &) const;
    bool AllowedAbsorption(Event &, size_t, size_t) const override;
};

} // namespace achilles

#endif
