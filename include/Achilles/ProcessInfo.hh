#ifndef PROCESS_INFO_HH
#define PROCESS_INFO_HH

#include "Achilles/ParticleInfo.hh"
#include "fmt/format.h"

#include <map>
#include <sstream>
#include <string>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

class PID;

struct ProcessInfo {
    using leptonic_state = std::pair<achilles::PID, std::vector<achilles::PID>>;
    using hadronic_state = std::pair<std::vector<achilles::PID>, std::vector<achilles::PID>>;
    using spectator_state = std::vector<achilles::PID>;
    ProcessInfo() = default;
    ProcessInfo(leptonic_state leptons) : m_leptonic{leptons} {}

    leptonic_state m_leptonic{};
    hadronic_state m_hadronic{};
    spectator_state m_spectator{};
    size_t FinalStateMultiplicity() const;
    size_t Multiplicity() const;
    std::vector<double> Masses() const;
    std::vector<long> Ids() const;
    std::vector<long> m_mom_map;
    int LeptonicCharge() const;

    bool operator==(const ProcessInfo &other) const {
        return m_leptonic == other.m_leptonic && m_hadronic == other.m_hadronic &&
               m_spectator == other.m_spectator;
        ;
    }

    template <typename OStream> friend OStream &operator<<(OStream &os, const ProcessInfo &info) {
        os << fmt::format("{}", info);
        return os;
    }

    bool SaveState(std::ostream &) const;
    bool LoadState(std::istream &);
};

} // namespace achilles

namespace fmt {

template <> struct formatter<achilles::ProcessInfo> {
    constexpr auto parse(format_parse_context &ctx) -> format_parse_context::iterator {
        return ctx.begin();
    }

    auto format(const achilles::ProcessInfo &info, format_context &ctx) const
        -> format_context::iterator {
        return format_to(
            ctx.out(), "ProcessInfo([{}, {}] -> [{}, {}])", info.m_leptonic.first,
            fmt::join(info.m_hadronic.first.begin(), info.m_hadronic.first.end(), ", "),
            fmt::join(info.m_leptonic.second.begin(), info.m_leptonic.second.end(), ", "),
            fmt::join(info.m_hadronic.second.begin(), info.m_hadronic.second.end(), ", "));
    }
};
} // namespace fmt

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

template <> struct std::hash<achilles::ProcessInfo> {
    std::size_t operator()(const achilles::ProcessInfo &p) const {
        return std::hash<std::string>{}(fmt::format("{}", p));
    }
};

#endif
