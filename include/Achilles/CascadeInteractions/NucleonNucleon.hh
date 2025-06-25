#ifndef ACHILLES_CASCADEINTERACTIONS_NUCLEONNUCLEON
#define ACHILLES_CASCADEINTERACTIONS_NUCLEONNUCLEON

#include "Achilles/DecayHandler.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/Units.hh"

namespace achilles {

class NucleonNucleon : public Interaction, RegistrableInteraction<NucleonNucleon> {
  private:
    enum class Mode {
        GiBUU = 1,
    };
    friend struct YAML::convert<NucleonNucleon::Mode>;
    enum class ResonanceMode {
        Decay = 1,
        Propagate = 2,
    };
    friend struct YAML::convert<NucleonNucleon::ResonanceMode>;

  public:
    NucleonNucleon();
    NucleonNucleon(const YAML::Node &);
    ~NucleonNucleon() = default;

    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<NucleonNucleon>(data);
    }

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return NucleonNucleon::Name(); }
    static std::string Name() { return "NucleonNucleon"; }

    // These functions are defined in the base class
    std::vector<std::pair<PID, PID>> InitialStates() const override;
    InteractionResults CrossSection(Event &, size_t, size_t) const override;
    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

    std::pair<double, double> TestNNElastic(double sqrts) const;

  private:
    Mode mode = NucleonNucleon::Mode::GiBUU;
    ResonanceMode resonance_mode = NucleonNucleon::ResonanceMode::Propagate;
    double NNElastic(double, PID, PID) const;
    const std::map<std::pair<PID, PID>, std::vector<PID>> outgoing;
    std::vector<std::pair<PID, PID>> AllowedResonanceStates(const Particle &particle1,
                                                            const Particle &particle2) const;
    double exp_sup = 0.;

    // Functions to cache integrals to speed up code
    double sqrts_min = 1800 / 1_GeV, sqrts_max = 6000 / 1_GeV;
    double mass_min = 1000 / 1_GeV, mass_max = 4000 / 1_GeV;
    static constexpr size_t nsqrts = 101, nmass = 101;
    void InitializeInterpolators();
    double SigmaNN2NDelta(double sqrts, double pcm, PID delta_id) const;
    double SigmaNN2NDeltaInterp(double sqrts, double pcm, PID delta_id) const;
    double DSigmaDMInterp(bool iresonance, double sqrts, double mdeltda, PID delta_id) const;
    Interp1D dsigma_ndelta;
    Interp2D dsigma_dm_ndelta;
    Interp2D dsigma_res_ndelta;
    DecayHandler decay_handler;
};

} // namespace achilles

namespace YAML {
template <> struct convert<achilles::NucleonNucleon::Mode> {
    static bool decode(const Node &node, achilles::NucleonNucleon::Mode &mode) {
        if(node.as<std::string>() == "GiBUU")
            mode = achilles::NucleonNucleon::Mode::GiBUU;
        else
            return false;
        return true;
    }
};

template <> struct convert<achilles::NucleonNucleon::ResonanceMode> {
    static bool decode(const Node &node, achilles::NucleonNucleon::ResonanceMode &mode) {
        if(node.as<std::string>() == "Decay")
            mode = achilles::NucleonNucleon::ResonanceMode::Decay;
        else if(node.as<std::string>() == "Propagate")
            mode = achilles::NucleonNucleon::ResonanceMode::Propagate;
        else
            return false;
        return true;
    }
};

} // namespace YAML

#endif
