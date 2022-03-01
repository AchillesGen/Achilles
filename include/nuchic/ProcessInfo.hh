#ifndef PROCESS_INFO_HH
#define PROCESS_INFO_HH

#include "nuchic/ParticleInfo.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace nuchic {
    using nuclear_map = std::map<std::vector<nuchic::PID>, std::vector<nuchic::PID>>;
    using nuclear_pair = std::pair<const std::vector<nuchic::PID>, std::vector<nuchic::PID>>;
}

namespace fmt {

template<>
struct formatter<nuchic::nuclear_pair> {
    template<typename ParseContext>
    constexpr auto parse(ParseContext &ctx) {
        return ctx.begin();
    }

    template<typename FormatContext>
    auto format(const nuchic::nuclear_pair &npair, FormatContext &ctx) {
        return format_to(ctx.out(), "[{}] -> [{}]",
                         join(npair.first.begin(), npair.first.end(), ", "),
                         join(npair.second.begin(), npair.second.end(), ", "));
    }
};

template<>
struct formatter<nuchic::nuclear_map> {
    template<typename ParseContext>
    constexpr auto parse(ParseContext &ctx) {
        return ctx.begin();
    }

    template<typename FormatContext>
    auto format(const nuchic::nuclear_map &nmap, FormatContext &ctx) {
        return format_to(ctx.out(), "{{{}}}", join(nmap.begin(), nmap.end(), ", ")); 
    }
};

}

namespace nuchic {

class Beam;

struct Process_Info {
    std::string m_model{};
    std::vector<nuchic::PID> m_ids{};
    nuclear_map m_states{};
    Process_Info() = default;
    Process_Info(std::string model, std::vector<nuchic::PID> ids={}) 
        : m_model(std::move(model)), m_ids(std::move(ids)) {}
    size_t Multiplicity() const;
    std::vector<double> Masses() const;
    std::vector<long> Ids() const;
    std::map<size_t, long> m_mom_map;

    bool operator==(const Process_Info &other) const {
        return m_model == other.m_model && m_ids == other.m_ids && m_states == other.m_states;
    }

    template<typename OStream>
    friend OStream& operator<<(OStream &os, const Process_Info &info) {
        os << fmt::format("Process_Info({}, PIDs = [{}], ",
                          info.m_model, fmt::join(info.m_ids.begin(), info.m_ids.end(), ", "));
        os << fmt::format("States = [{}])", info.m_states);
        return os;
    }
};

}

namespace YAML {

template<>
struct convert<nuchic::Process_Info> {
    static bool decode(const Node &node, nuchic::Process_Info &info) {
        if(!node.IsMap()) return false;
        if(!node["Model"].IsScalar()) return false;
        std::string model = node["Model"].as<std::string>();
        if(!node["Final States"].IsSequence()) return false;
        info = nuchic::Process_Info(model, node["Final States"].as<std::vector<nuchic::PID>>());
        return true;
    }
};

}

#endif
