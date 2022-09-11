#ifndef PROCESS_INFO_HH
#define PROCESS_INFO_HH

#include "Achilles/ParticleInfo.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

using nuclear_state = std::pair<std::vector<achilles::PID>, std::vector<achilles::PID>>;
class Beam;

struct Process_Info {
    std::vector<achilles::PID> ids{};
    nuclear_state state{};
    std::map<size_t, long> m_mom_map{};

    Process_Info(std::vector<achilles::PID> _ids={}) 
        : ids(std::move(_ids)) {}
    size_t Multiplicity() const;
    std::vector<double> Masses() const;
    std::vector<long> Ids() const;

    bool operator==(const Process_Info &other) const {
        return ids == other.ids && state == other.state;
    }

    template<typename OStream>
    friend OStream& operator<<(OStream &os, const Process_Info &info) {
        os << fmt::format("Process_Info(PIDs = [{}], ",
                          fmt::join(info.ids.begin(), info.ids.end(), ", "));
        os << fmt::format("State = [{}])", info.state);
        return os;
    }
};

struct Process_Group {
    std::string model{};
    size_t Multiplicity() const;
    size_t nucleons_in{}, nucleons_out{}, leptons_out{};
    static constexpr size_t leptons_in = 1;
    void AddProcess(const Process_Info &info);
    const std::vector<Process_Info> &Processes() const { return processes; }
    Process_Info& Process(size_t idx) { return processes[idx]; }
    const Process_Info& Process(size_t idx) const { return processes[idx]; }
    Process_Info& back() { return processes.back(); }
    Process_Info& front() { return processes.front(); }

    template<typename OStream>
    friend OStream& operator<<(OStream &os, const Process_Group &group) {
        os << fmt::format("Process_Group(Model = {}, Processes = [{}])",
                          group.model,
                          fmt::join(group.processes.begin(), group.processes.end(), ", "));
        return os;
    }

    private:
        std::vector<Process_Info> processes;
};

}

namespace fmt {

template<>
struct formatter<achilles::nuclear_state> {
    template<typename ParseContext>
    constexpr auto parse(ParseContext &ctx) {
        return ctx.begin();
    }

    template<typename FormatContext>
    auto format(const achilles::nuclear_state &npair, FormatContext &ctx) {
        return format_to(ctx.out(), "[{}] -> [{}]",
                         join(npair.first.begin(), npair.first.end(), ", "),
                         join(npair.second.begin(), npair.second.end(), ", "));
    }
};

template<>
struct formatter<achilles::Process_Info> {
    template<typename ParseContext>
    constexpr auto parse(ParseContext &ctx) {
        return ctx.begin();
    }

    template<typename FormatContext>
    auto format(const achilles::Process_Info &info, FormatContext &ctx) {
        return format_to(ctx.out(), "Process_Info(PIDs = [{}], States = [{}])",
                         fmt::join(info.ids.begin(), info.ids.end(), ", "));
    }
};

template<>
struct formatter<achilles::Process_Group> {
    template<typename ParseContext>
    constexpr auto parse(ParseContext &ctx) {
        return ctx.begin();
    }

    template<typename FormatContext>
    auto format(const achilles::Process_Group &group, FormatContext &ctx) {
        return format_to(ctx.out(), "Process_Group(Model = {}, Processes = [{}])",
                         group.model, fmt::join(group.Processes().begin(), group.Processes().end(), ", "));
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
        info = achilles::Process_Info(node["Final States"].as<std::vector<achilles::PID>>());
        return true;
    }
};

}

#endif
