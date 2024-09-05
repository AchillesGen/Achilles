#ifndef ACHILLES_DECAY_HANDLER
#define ACHILLES_DECAY_HANDLER

#include "Achilles/ParticleInfo.hh"

#include <map>
#include <string>
#include <vector>

namespace achilles {

class Particle;

struct DecayMode {
    double branching_ratio;
    size_t angular_mom;
    std::vector<PID> out_ids;
};

class DecayHandler {
  public:
    DecayHandler() = default;
    DecayHandler(const std::string &, double tolerance = 1e-8);
    std::vector<Particle> Decay(const Particle &) const;
    std::vector<DecayMode> AllowedDecays(PID) const;
    double BranchingRatio(PID, std::vector<PID>) const;

  private:
    std::vector<double> BranchingRatios(const Particle &) const;
    std::vector<Particle> TwoBodyDecay(double m2, const std::vector<PID> &, size_t) const;
    std::map<PID, std::vector<DecayMode>> m_decays;
};

} // namespace achilles

#endif
