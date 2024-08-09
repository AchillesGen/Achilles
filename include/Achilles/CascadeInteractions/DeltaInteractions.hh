#ifndef ACHILLES_CASCADEINTERACTIONS_DELTAINTERACTIONS
#define ACHILLES_CASCADEINTERACTIONS_DELTAINTERACTIONS

#include "Achilles/Interactions.hh"

namespace achilles {

class DeltaInteraction : public Interaction, RegistrableInteraction<DeltaInteraction> {
  private:
    enum class Mode {
        GiBUU = 1,
    };
    friend struct YAML::convert<DeltaInteraction::Mode>;

  public:
    DeltaInteraction() = default;
    DeltaInteraction(const YAML::Node &);

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
    double TestInelastic1Pi(double sqrts, size_t type) const;
    double TestDeltaDSigma(bool iresonance, double sqrts, double mdelta) const;
    double TestDeltaDSigmaDOmegaDM(double cost, double sqrts, double mdelta, PID);
    double TestDeltaSigma(double sqrts) const;

  private:
    Mode mode;

    double NNElastic(double, PID, PID) const;

    // Calculation based on https://doi.org/10.1016/0375-9474(94)90405-7
    const std::map<std::pair<PID, PID>, double> sigma_max;
    const std::map<std::pair<PID, PID>, std::vector<PID>> outgoing;

    static constexpr double beta = 300;
    double vfunc(double q) const { return beta * beta / (beta * beta + q * q); }

    struct NNInelastic1PiParams {
        double norm, a, b, n1, n2;
    };
    static constexpr NNInelastic1PiParams pppi0{61.3, 1.52, 2.50, 6.18, 3.48};
    static constexpr NNInelastic1PiParams pnpip{122.6, 1.52, 2.50, 6.18, 3.48};
    static constexpr NNInelastic1PiParams pppim{24.9, 3.30, 0.85, 1.93, 0.002};
    static constexpr NNInelastic1PiParams pnpi0{7.25, 0.88, 0, 2.31, 3.64};
    static constexpr std::array<NNInelastic1PiParams, 4> inelastic_1pi_params{pppi0, pnpip, pppim,
                                                                              pnpi0};
    double NNInelastic1PiBackground(double, const NNInelastic1PiParams &) const;
    double NNResonance1Pi(double, PID, PID, PID) const;

    // NN to NDelta and NDelta to NN
    double pcm2(double, double, double) const;
    double NN2NDelta(const Particle &p1, const Particle &p2);
    double NDelta2NN(const Particle &p1, const Particle &p2);
    double SigmaNN2NDelta(double sqrts, PID delta_id) const;
    double DSigmaDM(bool iresonance, double sqrts, double mdelta, PID delta_id) const;
    double MatNN2NDelta(double t, double u, double mdelta, PID delta_id) const;

    // TODO: Should be moved out of the class
    double GetEffectiveWidth(PID id, double mass, double mass1, double mass2,
                             size_t angular_mom) const;
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
