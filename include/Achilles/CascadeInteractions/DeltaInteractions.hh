#ifndef ACHILLES_CASCADEINTERACTIONS_DELTAINTERACTIONS
#define ACHILLES_CASCADEINTERACTIONS_DELTAINTERACTIONS

#include "Achilles/DecayHandler.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/Units.hh"

namespace achilles {

class DeltaInteraction : public Interaction, RegistrableInteraction<DeltaInteraction> {
  private:
    using absorption_states = std::vector<std::pair<size_t, std::vector<PID>>>;
    enum class Mode {
        GiBUU = 1,
    };
    friend struct YAML::convert<DeltaInteraction::Mode>;

  public:
    DeltaInteraction();
    DeltaInteraction(const YAML::Node &);
    ~DeltaInteraction();

    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<DeltaInteraction>(data);
    }

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return DeltaInteraction::Name(); }
    static std::string Name() { return "DeltaInteraction"; }

    // These functions are defined in the base class
    std::vector<std::pair<PID, PID>> InitialStates() const override;
    InteractionResults CrossSection(Event &, size_t, size_t) const override;
    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

    std::pair<double, double> TestNNElastic(double sqrts) const;
    double TestDeltaDSigma(bool iresonance, double sqrts, double mdelta) const;
    double TestDeltaDSigmaDOmegaDM(double cost, double sqrts, double mdelta, PID);
    double TestDeltaSigma(double sqrts) const;
    double TestNPiSigma(const Particle &p1, const Particle &p2, PID res) const;
    double TestNDelta2NDelta(const Particle &p1, const Particle &p2, PID dout, PID nout) const;

  private:
    Mode mode;
    const std::map<std::pair<PID, PID>, std::vector<PID>> outgoing;

    // NDelta to NN
    double SigmaNDelta2NN(double sqrts, double pcm, PID delta_id, PID nucleon, double mass) const;

    // Npi to X
    double SigmaNPi2Delta(const Particle &, const Particle &, PID) const;
    double SigmaNPi2N(const Particle &, const Particle &, PID) const;

    // NDelta to NDelta
    double SigmaNDelta2NDelta(const Particle &, const Particle &, PID, PID) const;
    double GetIso(int, int, int, int) const;

    // S-wave absorption cross section
    bool swave_enabled = true;
    absorption_states AllowedAbsorption(Event &, size_t, size_t) const;
    struct AbsorptionStates {
        PID absorption_partner;
        std::vector<PID> outgoing;
    };
    using AbsorptionModes = std::map<std::pair<PID, PID>, std::vector<AbsorptionStates>>;
    static const AbsorptionModes absorption_modes;
    mutable std::vector<Particle> absorption_partners;
    std::pair<size_t, size_t> FindClosest(Event &, size_t, size_t) const;
    std::vector<Particle> HandleAbsorption(const Particle &particle1, const Particle &particle2,
                                           const std::vector<PID> &out_pids, Random &ran) const;
    double exp_sup = 0.;

    // Functions to cache integrals to speed up code
    double sqrts_min = 1800 / 1_GeV, sqrts_max = 6000 / 1_GeV;
    double mass_min = 1000 / 1_GeV, mass_max = 4000 / 1_GeV;
    static constexpr size_t nsqrts = 101, nmass = 101;
    double DSigmaDMInterp(bool iresonance, double sqrts, double mdeltda, PID delta_id) const;
    Interp2D dsigma_dm_ndelta;
    Interp2D dsigma_res_ndelta;
    DecayHandler decay_handler;
};

} // namespace achilles

namespace YAML {
template <> struct convert<achilles::DeltaInteraction::Mode> {
    static bool decode(const Node &node, achilles::DeltaInteraction::Mode &mode) {
        if(node.as<std::string>() == "GiBUU") mode = achilles::DeltaInteraction::Mode::GiBUU;
        return true;
    }
};

} // namespace YAML

#endif
