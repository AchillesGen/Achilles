#include "Achilles/System.hh"
#include "Achilles/Utilities.hh"

std::vector<std::string> achilles::Filesystem::GetPluginPaths() {
    std::vector<std::string> dirs;
    const char *env = std::getenv("ACHILLES_PLUGIN_PATH");
    if(env) achilles::tokenize(std::string(env), dirs, ":");
    dirs.push_back(achilles::PathVariables::installLibs);
    return dirs;
}
