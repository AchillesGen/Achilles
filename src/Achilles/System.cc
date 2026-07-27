// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/System.hh"
#include "Achilles/Exception.hh"
#include "Achilles/Utilities.hh"
#include "spdlog/spdlog.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#include <regex>
#pragma GCC diagnostic pop

using achilles::Filesystem::Cache;

bool achilles::SystemVariables::IsPluginName(const std::string &filename) {
    static const std::regex validPluginName{std::string(pluginBase) + ".*\\" +
                                            std::string(libSuffix)};
    return std::regex_search(filename, validPluginName);
}

achilles::PathVariables &achilles::PathVariables::Instance() {
    static PathVariables global_paths;
    return global_paths;
}

// Function to load the search path for Achilles files
std::vector<fs::path> achilles::Filesystem::AchillesPath() {
    std::vector<fs::path> dirs;
    dirs.push_back(fs::current_path());
    const char *env = std::getenv("ACHILLES_PATH");
    std::vector<std::string> tmp;
    if(env) tokenize(std::string(env), tmp, ":");
    for(const auto &path : tmp) dirs.push_back(fs::path(path));
    // Shared data cache populated by the Python downloader. The directory name is
    // chosen entirely on the Python side (cross-platform) and handed to us via the
    // environment. Falls back to SharePath() for standalone installs.
    const char *dataDir = std::getenv("ACHILLES_DATA_DIR");
    if(dataDir) dirs.push_back(fs::path(dataDir));
    dirs.push_back(PathVariables::Instance().SharePath());
    return dirs;
}

std::vector<std::string> achilles::Filesystem::GetPluginPaths() {
    std::vector<std::string> dirs;
    const char *env = std::getenv("ACHILLES_PLUGIN_PATH");
    if(env) achilles::tokenize(std::string(env), dirs, ":");
    dirs.push_back(achilles::PathVariables::Instance().LibsPath());
    return dirs;
}

std::string achilles::Filesystem::FindFile(const std::string &filename, const std::string &head) {
    spdlog::trace("{}: Loading file {}", head, filename);
    auto dirs = AchillesPath();

    for(const auto &path : dirs) {
        if(fs::exists(path / filename)) {
            spdlog::debug("{}: Found {} at {}", head, filename, path);
            return path / filename;
        }
        spdlog::debug("{}: Could not find {} at {}", head, filename, path);
    }

    spdlog::warn("{}: Could not load {} from {}", head, filename,
                 fmt::join(dirs.begin(), dirs.end(), ":"));
    throw AchillesLoadError(filename);
}

std::string achilles::Filesystem::FindFlux(const std::string &filename, const std::string &head) {
    spdlog::trace("{}: Loading flux file {}", head, filename);
    auto dirs = AchillesPath();
    dirs.push_back(PathVariables::Instance().FluxPath());

    for(const auto &path : dirs) {
        if(fs::exists(path / filename)) {
            spdlog::debug("{}: Found {} at {}", head, filename, path);
            return path / filename;
        }
        spdlog::debug("{}: Could not find {} at {}", head, filename, path);
    }

    spdlog::warn("{}: Could not load {} from {}", head, filename,
                 fmt::join(dirs.begin(), dirs.end(), ":"));
    throw AchillesLoadError(filename);
}

bool Cache::FindCachedState(std::size_t hash) {
    auto cachePath = Path() / fmt::format("{:x}", hash);
    return fs::exists(cachePath);
}
