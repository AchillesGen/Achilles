#include "spdlog/spdlog.h"
#include <algorithm>

#include "Achilles/Interactions.hh"

#include "plugins/InteractionLoader.hh"
#include "plugins/System.hh"

#include <dlfcn.h>

using achilles::InteractionLoader;
using namespace achilles::SystemVariables;

bool InteractionLoader::m_init{false};
std::vector<fs::path> InteractionLoader::plugins{};
std::set<achilles::Interactions *> loadPlugins;

bool InteractionLoader::LoadInteractions(const std::vector<std::string> &paths) {
    // Ensure that the plugins are only loaded once
    if(m_init) return true;

    // Find all plugins in the paths provided
    for(const auto &path : paths) {
        spdlog::info("Loading plugins found at: {}", path);
        if(fs::exists(path) && fs::is_directory(path)) {
            // Loop over all files in the directory
            for(const auto &entry : fs::directory_iterator(path)) {
                auto filename = entry.path().filename().string();
                if(filename.find(librarySuffix) == std::string::npos) continue;
                if(std::find(plugins.begin(), plugins.end(), filename) == plugins.end())
                    plugins.push_back(entry);
            }
        }
    }

    // Load plugins
    for(const auto &plugin : plugins) {
        spdlog::info("Loading {}...", plugin.string());
        void *handle = dlopen(plugin.c_str(), RTLD_NOW);
        if(!handle) { spdlog::warn("Cannot open {}: {}", plugin, dlerror()); }
    }

    return true;
}
