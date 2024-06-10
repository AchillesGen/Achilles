#include "plugins/Manager/PluginManager.hh"
#include "Achilles/System.hh"
#include "Achilles/Version.hh"
#include "fmt/ranges.h"
#include "spdlog/spdlog.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#include <regex>
#pragma GCC diagnostic pop
#include <dlfcn.h>

using achilles::Plugin::Manager;

template <typename Ret, typename... Args> std::function<Ret(Args...)> fptr_cast(void *fptr) {
    using function_type = Ret (*)(Args...);
    using ptr_size_type = std::conditional_t<sizeof(fptr) == 4, long, long long>;

    return reinterpret_cast<function_type>(reinterpret_cast<ptr_size_type>(fptr));
}

void Manager::Open() {
    auto paths = achilles::Filesystem::GetPluginPaths();
    const auto load_plugin = [&](const auto &file) {
        auto filename = file.path().filename().string();
        if(std::regex_search(filename, SystemVariables::validPluginName)) {
            void *handle = dlopen(file.path().c_str(), RTLD_NOW);
            spdlog::debug("PluginManager: Loading plugin {}", filename);
            if(!handle) {
                spdlog::warn("PluginManager: Could not load plugin: {}", dlerror());
                return;
            }
            auto get_version =
                fptr_cast<void, std::array<int, 3> &>(dlsym(handle, "ExpectedVersion"));

            if(dlerror() != NULL) {
                spdlog::warn("PluginManager: Plugin does not define an Achilles version");
                dlclose(handle);
                return;
            }

            std::array<int, 3> expected_version;
            get_version(expected_version);
            if(!ValidateVersion(expected_version)) {
                spdlog::warn("PluginManager: Plugin expects version {}, using version {}",
                             fmt::join(expected_version.begin(), expected_version.end(), "."),
                             fmt::join(CurrentVersion.begin(), CurrentVersion.end(), "."));
                dlclose(handle);
                return;
            }

            auto register_plugin = fptr_cast<void>(dlsym(handle, "Register"));
            if(register_plugin) register_plugin();

            m_handles.push_back(std::move(handle));
        }
    };

    for(const auto &dir : paths) {
        const auto &directory{dir};
        Filesystem::RecurseDirectories(directory, load_plugin,
                                       achilles::Filesystem::FileType::regular);
    }
}

void Manager::Close() {
    for(auto &handle : m_handles) {
        if(handle) dlclose(handle);
    }
}

bool Manager::ValidateVersion(const std::array<int, 3> &version) const {
    return version[0] == CurrentVersion[0] && version[1] == CurrentVersion[1];
}
