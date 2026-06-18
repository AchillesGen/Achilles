#ifndef ACHILLES_CASCADEINTERACTIONS_NASAINTERACTIONS
#define ACHILLES_CASCADEINTERACTIONS_NASAINTERACTIONS

#include "Achilles/Interactions.hh"

namespace achilles {

class ThreeVector;

class NasaInteraction : public Interaction, RegistrableInteraction<NasaInteraction> {
  public:
    NasaInteraction() = default;
    NasaInteraction(const YAML::Node &) {};

    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<NasaInteraction>(data);
    }

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return NasaInteraction::Name(); }
    static std::string Name() { return "NasaInteraction"; }

    // These functions are defined in the base class
    std::vector<std::pair<PID, PID>> InitialStates() const override;
    InteractionResults CrossSection(Event &, size_t, size_t) const override;
    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

  private:
    ThreeVector MakeMomentum(bool, double, const std::vector<double> &) const;
};

} // namespace achilles

#endif
