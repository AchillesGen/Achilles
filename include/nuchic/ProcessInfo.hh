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
    std::vector<nuchic::PID> m_ids;
    Process_Info(std::string model, std::vector<nuchic::PID> ids={}) 
        : m_model(std::move(model)), m_ids(std::move(ids)) {}
    // void AddBeam(const nuchic::Beam&);

    template<typename OStream>
    friend OStream& operator<<(OStream &os, const Process_Info &info) {
        os << "Process_Info(" << info.m_model;
        for(const auto &pid : info.m_ids) {
            os << ", " << static_cast<int>(pid);
        }
        os << ")";
        return os;
    }
};

}

namespace YAML {

template<>
struct convert<std::vector<nuchic::Process_Info>> {
    static bool decode(const Node &node, std::vector<nuchic::Process_Info> &info) {
        if(!node.IsMap()) return false;
        if(!node["Model"].IsScalar()) return false;
        std::string model = node["Model"].as<std::string>();
        if(!node["Leptons"].IsSequence()) return false;
        for(auto subnode : node["Leptons"]) {
            if(!subnode.IsSequence()) return false;
            info.emplace_back(model,
                              subnode.as<std::vector<nuchic::PID>>());
        }
        return true;
    }
};

}

#endif
