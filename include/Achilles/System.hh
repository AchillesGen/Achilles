#ifndef SYSTEM_HH
#define SYSTEM_HH

#include <string>

namespace achilles {

namespace SystemVariables {
    const std::string libPrefix = "lib";
    const std::string libSuffix = ".so";
}

namespace PathVariables {
    const std::string buildPath = "/home/isaacson/Documents/Projects/NeutrinoGenerator/Achilles/optional_bsm/build/";
    const std::string installPath = "/usr/local/";
    const std::string buildLibs = buildPath + "lib/";
    const std::string installLibs = installPath + "lib/";
    const std::string buildData = buildPath + "data/";
    const std::string installData = installPath + "share/achilles/data/";
}

}

#endif
