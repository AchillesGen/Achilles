#ifndef RUNMODES_HH
#define RUNMODES_HH

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace nuchic {

enum class RunMode {
    FullPhaseSpace,
    FixedAngle,
    FixedAngleEnergy
};

}

namespace YAML {

template<>
struct convert<nuchic::RunMode> {
    static bool decode(const Node &node, nuchic::RunMode &rhs) {
        if(node.as<std::string>() == "FullPhaseSpace") {
            rhs = nuchic::RunMode::FullPhaseSpace;
            return true;
        } else if(node.as<std::string>() == "FixedAngle") {
            rhs = nuchic::RunMode::FixedAngle;
            return true;
        } else if(node.as<std::string>() == "FixedAngleEnergy") {
            rhs = nuchic::RunMode::FixedAngleEnergy;
            return true;
        }

        return false;
    }
};

}

#endif
