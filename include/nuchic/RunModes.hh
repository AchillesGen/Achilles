#ifndef RUNMODES_HH
#define RUNMODES_HH

#include "yaml-cpp/yaml.h"

namespace nuchic {

enum class RunMode {
    FullPhaseSpace,
    FixedAngle
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
        }

        return false;
    }
};

}

#endif
