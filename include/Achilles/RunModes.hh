#ifndef RUNMODES_HH
#define RUNMODES_HH

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

enum class RunMode { FullPhaseSpace, FixedAngle, FixedAngleEnergy };

}

namespace YAML {

template <> struct convert<achilles::RunMode> {
    static bool decode(const Node &node, achilles::RunMode &rhs) {
        if(node.as<std::string>() == "FullPhaseSpace") {
            rhs = achilles::RunMode::FullPhaseSpace;
            return true;
        } else if(node.as<std::string>() == "FixedAngle") {
            rhs = achilles::RunMode::FixedAngle;
            return true;
        } else if(node.as<std::string>() == "FixedAngleEnergy") {
            rhs = achilles::RunMode::FixedAngleEnergy;
            return true;
        }

        return false;
    }
};

} // namespace YAML

#endif
