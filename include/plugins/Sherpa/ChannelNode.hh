#ifndef PLUGINS_SHERPA_CHANNELNODE
#define PLUGINS_SHERPA_CHANNELNODE

#include "Achilles/ParticleInfo.hh"

#include <map>
#include <set>

namespace achilles::channels {

struct ParticleInfo {
    unsigned int idx = 0;
    achilles::ParticleInfo info{0};
};

using Current = std::unordered_map<unsigned int, std::set<ParticleInfo>>;
using DecayProds = std::pair<ParticleInfo, ParticleInfo>;
using DecayMap = std::map<ParticleInfo, DecayProds>;
using DecayChain = std::map<ParticleInfo, std::set<DecayProds>, std::greater<ParticleInfo>>;

struct DataFrame {
    std::set<int> pid;
    std::vector<unsigned int> avail_currents;
    std::vector<unsigned int> currents{};
    unsigned int idx_sum{};
};

struct Description {
    std::vector<ParticleInfo> info;
    DecayMap decays{};
};
using ChannelVec = std::vector<Description>;

ChannelVec AddDecays(unsigned int, const Description &, const DecayChain &);
bool operator<(const ParticleInfo &a, const ParticleInfo &b);
bool operator==(const ParticleInfo &a, const ParticleInfo &b);
inline bool operator!=(const ParticleInfo &a, const ParticleInfo &b) {
    return !(a == b);
}
inline bool operator>(const ParticleInfo &a, const ParticleInfo &b) {
    return b < a && b != a;
}
bool operator<(const Description &a, const Description &b);

inline std::string ToString(ParticleInfo info) {
    return fmt::format("{}", info.info.Name());
}

inline std::string ToString(Description desc) {
    std::vector<std::string> particles;
    for(const auto &info : desc.info) particles.push_back(ToString(info));
    std::vector<std::string> decays;
    for(const auto &decay : desc.decays) {
        decays.push_back(fmt::format("Decay({} -> ({}, {}))", ToString(decay.first),
                                     ToString(decay.second.first), ToString(decay.second.second)));
    }
    return fmt::format("ChannelDescription({}, {})", fmt::join(particles, ", "),
                       fmt::join(decays, ", "));
}

std::string ToString(DecayChain chain);

} // namespace achilles::channels

#endif // PLUGINS_SHERPA_CHANNELNODE
