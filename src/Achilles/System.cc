#include "Achilles/Exception.hh"
#include "Achilles/System.hh"
#include "Achilles/Utilities.hh"
#include "spdlog/spdlog.h"
#include "yaml-cpp/yaml.h"
#include "fmt/std.h"

// Function to load the search path for Achilles files
std::vector<fs::path> achilles::Filesystem::AchillesPath() {
    std::vector<fs::path> dirs; 
    dirs.push_back(fs::current_path());
    const char *env = std::getenv("ACHILLES_PATH");
    std::vector<std::string> tmp;
    if(env) tokenize(std::string(env), tmp, ":");
    for(const auto &path : tmp) dirs.push_back(fs::path(path));
    dirs.push_back(PathVariables::installShare);
    return dirs;
}

std::vector<std::string> achilles::Filesystem::GetPluginPaths() {
    std::vector<std::string> dirs;
    const char* env = std::getenv("ACHILLES_PLUGIN_PATH");
    if(env) achilles::tokenize(std::string(env), dirs, ":");
    dirs.push_back(achilles::PathVariables::installLibs);
    return dirs;
}

std::string achilles::Filesystem::FindFile(const std::string &filename, const std::string &head) {
    auto dirs = AchillesPath();

    for(const auto &path : dirs) {
        if(fs::exists(path / filename)) {
            spdlog::debug("{}: Found {} at {}", head, filename, path);
            return path / filename;
        }
        spdlog::debug("{}: Could not find {} at {}", head, filename, path);
    }

    spdlog::warn("{}: Could not load {} from {}", head, filename, fmt::join(dirs.begin(), dirs.end(), ":"));
    throw AchillesLoadError(filename);
}
