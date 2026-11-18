// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/System.hh"
#include "Achilles/Exception.hh"
#include "Achilles/Utilities.hh"
#include "spdlog/spdlog.h"

#include <dlfcn.h>

using achilles::Filesystem::Cache;

const std::regex &achilles::SystemVariables::ValidPluginName() {
    static const std::regex valid(std::string(ACHILLES_LIB_PREFIX) + "AchillesPlugin.*\\" +
                                  ACHILLES_LIB_SUFFIX);
    return valid;
}

achilles::PathVariables &achilles::PathVariables::Instance() {
    static PathVariables global_paths;
    return global_paths;
}

void achilles::PathVariables::InitializeStandalone() {
    // Any address inside this library identifies the file it was loaded from
    static const char anchor{};
    Dl_info info{};
    if(dladdr(&anchor, &info) == 0 || info.dli_fname == nullptr) return;

    std::error_code ec;
    auto library = fs::canonical(fs::path(info.dli_fname), ec);
    if(ec) return;

    libsPath = library.parent_path();
    sharePath = (libsPath / ACHILLES_LIB_TO_SHARE).lexically_normal();
    rootPath = (libsPath / ACHILLES_LIB_TO_ROOT).lexically_normal();
    fluxPath = sharePath / "flux";
    dataPath = sharePath / "data";
    defaultsPath = dataPath / "default";
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
