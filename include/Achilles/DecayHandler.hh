#ifndef ACHILLES_DECAY_HANDLER
#define ACHILLES_DECAY_HANDLER

#include "Achilles/ParticleInfo.hh"

#include <fstream>
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

    // TEST: angular distribution
    static std::ofstream m_angular_dist_file;
};

} // namespace achilles

namespace fmt {

template <> struct formatter<achilles::DecayMode> {
    template <typename ParseContext> constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const achilles::DecayMode &mode, FormatContext &ctx) const {
        return format_to(ctx.out(), "Mode[{{{}}}, Br={}]", fmt::join(mode.out_ids, ", "),
                         mode.branching_ratio);
    }
};

} // namespace fmt

#endif
