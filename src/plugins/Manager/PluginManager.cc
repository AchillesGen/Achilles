#include "plugins/Manager/PluginManager.hh"
#include "spdlog/spdlog.h"
#include <dlfcn.h>
#include <filesystem>

using achilles::Plugin::Manager;

std::vector<std::string> Manager::GetPluginPaths() {
    std::vector<std::string> dirs;
    const char* env = std::getenv("ACHILLES_PLUGIN_PATH");
    if(env) achilles::tokenize(std::string(env), dirs, ":");
    dirs.push_back(achilles::PathVariables::installLibs);
    return dirs;
}

void Manager::Open() {
    auto paths = GetPluginPaths();
    for(const auto &dir : paths) {
        std::filesystem::path directory{dir};
        for(const auto &file_path : std::filesystem::recursive_directory_iterator(directory)) {
            if(file_path.is_regular_file()) {
                auto file = file_path.path().string();
                if(file.find("AchillesPlugin") != std::string::npos
                        && file.find(achilles::SystemVariables::libPrefix) != std::string::npos
                        && file.find(achilles::SystemVariables::libSuffix) != std::string::npos) {
                    void *handle = dlopen(file.c_str(), RTLD_NOW);
                    m_handles.push_back(handle);
                    spdlog::info("PluginManager: Loading plugin {}", file);
                    if(!handle) {
                        spdlog::warn("PluginManager: Could not load plugin: {}", dlerror());
                    }
                }
            }
        }
    }
}

void Manager::Close() {
    for(auto &handle : m_handles) {
        if(handle) dlclose(handle);
    }
}
