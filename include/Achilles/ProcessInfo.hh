#ifndef PROCESS_INFO_HH
#define PROCESS_INFO_HH

#include "Achilles/ParticleInfo.hh"
#include "fmt/format.h"

#include <map>
#include <sstream>
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
    std::map<size_t, long> m_mom_map;
    int LeptonicCharge() const;

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

namespace fmt {

template<>
struct formatter<achilles::Process_Info> {
    constexpr auto parse(format_parse_context &ctx) -> format_parse_context::iterator {
        return ctx.begin();
    }

    auto format(const achilles::Process_Info &info, format_context &ctx) const -> format_context::iterator {
        return format_to(ctx.out(), "Process_Info({}, PIDs = [{}], States = [{}])",
                         info.m_model, fmt::join(info.m_ids.begin(), info.m_ids.end(), ", "), info.m_states);
    }
};
}

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

namespace fmt {

template <> struct formatter<achilles::ProcessInfo> {
    template <typename ParseContext> constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const achilles::ProcessInfo &procinfo, FormatContext &ctx) const {
        std::stringstream ss;
        ss << procinfo;

        return format_to(ctx.out(), ss.str());
    }
};

} // namespace fmt

#endif
