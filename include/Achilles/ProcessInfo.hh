#ifndef PROCESS_INFO_HH
#define PROCESS_INFO_HH

#include "Achilles/ParticleInfo.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {
    using nuclear_map = std::map<std::vector<achilles::PID>, std::vector<achilles::PID>>;
    using nuclear_pair = std::pair<const std::vector<achilles::PID>, std::vector<achilles::PID>>;
}

namespace fmt {

template<>
struct formatter<achilles::nuclear_pair> {
    constexpr auto parse(format_parse_context &ctx) -> format_parse_context::iterator {
        return ctx.begin();
    }

    auto format(const achilles::nuclear_pair &npair, format_context &ctx) const -> format_context::iterator {
        return format_to(ctx.out(), "[{}] -> [{}]",
                         join(npair.first.begin(), npair.first.end(), ", "),
                         join(npair.second.begin(), npair.second.end(), ", "));
    }
};

template<>
struct formatter<achilles::nuclear_map> {
    constexpr auto parse(format_parse_context &ctx) -> format_parse_context::iterator {
        return ctx.begin();
    }

    auto format(const achilles::nuclear_map &nmap, format_context &ctx) const -> format_context::iterator {
        return format_to(ctx.out(), "{{{}}}", join(nmap.begin(), nmap.end(), ", ")); 
    }
};

}

namespace achilles {

class Beam;

struct Process_Info {
    std::string m_model{};
    std::vector<achilles::PID> m_ids{};
    nuclear_map m_states{};
    Process_Info() = default;
    Process_Info(std::string model, std::vector<achilles::PID> ids={}) 
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

template<>
struct convert<achilles::Process_Info> {
    static bool decode(const Node &node, achilles::Process_Info &info) {
        if(!node.IsMap()) return false;
        if(!node["Model"].IsScalar()) return false;
        std::string model = node["Model"].as<std::string>();
        if(!node["Final States"].IsSequence()) return false;
        info = achilles::Process_Info(model, node["Final States"].as<std::vector<achilles::PID>>());
        return true;
    }
};

}

#endif
