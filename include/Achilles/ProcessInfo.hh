#ifndef PROCESS_INFO_HH
#define PROCESS_INFO_HH

#include "Achilles/ParticleInfo.hh"
#include "fmt/format.h"

#include <map>
#include <string>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

class PID;

struct ProcessInfo {
    using leptonic_state = std::pair<achilles::PID, std::vector<achilles::PID>>;
    using hadronic_state = std::pair<std::vector<achilles::PID>, std::vector<achilles::PID>>;
    ProcessInfo() = default;
    ProcessInfo(leptonic_state leptons) : m_leptonic{leptons} {}

    leptonic_state m_leptonic{};
    hadronic_state m_hadronic{};
    size_t Multiplicity() const;
    std::vector<double> Masses() const;
    std::vector<long> Ids() const;
    std::map<size_t, long> m_mom_map;
    int LeptonicCharge() const;
    size_t NInitialStates(size_t, size_t) const;

    bool operator==(const ProcessInfo &other) const {
        return m_leptonic == other.m_leptonic && m_hadronic == other.m_hadronic;
    }

    template <typename OStream> friend OStream &operator<<(OStream &os, const ProcessInfo &info) {
        os << fmt::format(
            "Process_Info([{}, {}] -> [{}, {}])", info.m_leptonic.first,
            fmt::join(info.m_hadronic.first.begin(), info.m_hadronic.first.end(), ", "),
            fmt::join(info.m_leptonic.second.begin(), info.m_leptonic.second.end(), ", "),
            fmt::join(info.m_hadronic.second.begin(), info.m_hadronic.second.end(), ", "));
        return os;
    }
};

} // namespace achilles

namespace YAML {

template <> struct convert<achilles::ProcessInfo> {
    static bool decode(const Node &node, achilles::ProcessInfo &info) {
        if(!node.IsMap()) return false;
        if(!node["Leptons"].IsSequence()) return false;
        const auto initial = node["Leptons"][0].as<achilles::PID>();
        const auto final = node["Leptons"][1].as<std::vector<achilles::PID>>();
        info = achilles::ProcessInfo({initial, final});
        return true;
    }
};

} // namespace YAML

#endif
