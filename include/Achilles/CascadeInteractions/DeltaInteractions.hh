#ifndef ACHILLES_CASCADEINTERACTIONS_DELTAINTERACTIONS
#define ACHILLES_CASCADEINTERACTIONS_DELTAINTERACTIONS

#include "Achilles/DecayHandler.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/Units.hh"

namespace achilles {

class DeltaInteraction : public Interaction, RegistrableInteraction<DeltaInteraction> {
  private:
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
    void TestInterpolation() const;

  private:
    Mode mode;

    double NNElastic(double, PID, PID) const;

    void printmode() const { spdlog::debug("mode = {}", static_cast<int>(mode)); }

    // Calculation based on https://doi.org/10.1016/0375-9474(94)90405-7
    const std::map<std::pair<PID, PID>, double> sigma_max;
    const std::map<std::pair<PID, PID>, std::vector<PID>> outgoing;

    static constexpr double beta = 300;
    double vfunc(double q) const { return beta * beta / (beta * beta + q * q); }

    // Npi to Delta
    double SigmaNPi2Delta(const Particle &, const Particle &, PID) const;

    // NN to NDelta and NDelta to NN
    double Pcm2(double, double, double) const;
    double SigmaNN2NDelta(double sqrts, double pcm, PID delta_id) const;
    double SigmaNDelta2NN(double sqrts, double pcm, PID delta_id, PID nucleon, double mdelta) const;
    double DSigmaDM(bool iresonance, double sqrts, double mdelta, PID delta_id) const;
    double MatNN2NDelta(double t, double u, double mdelta, PID delta_id) const;

    // NDelta to NDelta
    double SigmaNDelta2NDelta(const Particle &, const Particle &, PID, PID) const;
    double GetIso(int, int, int, int) const;
    double SpectralDelta(PID, double mu2) const;
    double DSigmaND2ND(double sqrts, double mn1, double mn2, double mu1, double mu2,
                       double spectral) const;
    double MatNDelta2NDelta(double t, double mu1, double mu2) const;

    // Functions to cache integrals to speed up code
    double sqrts_min = 1800 / 1_GeV, sqrts_max = 6000 / 1_GeV;
    double mass_min = 1000 / 1_GeV, mass_max = 4000 / 1_GeV;
    static constexpr size_t nsqrts = 101, nmass = 101;
    void InitializeInterpolators();
    double SigmaNN2NDeltaInterp(double sqrts, double pcm, PID delta_id) const;
    double DSigmaDMInterp(bool iresonance, double sqrts, double mdeltda, PID delta_id) const;
    Interp1D dsigma_ndelta, dsigma_ndnd;
    Interp2D dsigma_dm_ndelta;
    Interp2D dsigma_res_ndelta;

    // TODO: Should be moved out of the class
    double GetEffectiveWidth(PID id, double mass, double mass1, double mass2,
                             size_t angular_mom) const;
    double GenerateMass(const Particle &p1, const Particle &p2, PID res, PID other, Random &ran,
                        double smax) const;
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
