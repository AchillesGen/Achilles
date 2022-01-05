#ifndef PROCESS_INFO_HH
#define PROCESS_INFO_HH

#include "nuchic/ParticleInfo.hh"
#include <iostream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace nuchic {

class Beam;

struct Process_Info {
    std::string m_model{};
    std::vector<nuchic::PID> m_ids{};
    std::map<std::vector<nuchic::PID>, std::vector<nuchic::PID>> m_states{};
    Process_Info() = default;
    Process_Info(std::string model, std::vector<nuchic::PID> ids={}) 
        : m_model(std::move(model)), m_ids(std::move(ids)) {}
    size_t Multiplicity() const;
    std::vector<double> Masses() const;
    std::vector<long> Ids() const;
    std::map<size_t, long> m_mom_map;

    template<typename OStream>
    friend OStream& operator<<(OStream &os, const Process_Info &info) {
        os << "Process_Info(" << info.m_model;
        os << ", PIDs = [";
        for(const auto &pid : info.m_ids)
            os << static_cast<int>(pid) << ", ";
        os << "\b\b], States = [";
        for(const auto &state : info.m_states) {
            os << "{[";
            for(const auto &initial : state.first)
                os << static_cast<int>(initial) << ", ";
            os << "\b\b] -> [";
            for(const auto &final : state.second)
                os << static_cast<int>(final) << ", ";
            os << "\b\b]}, ";
        }
        os << "\b\b)";
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
