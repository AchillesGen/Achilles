#include "Achilles/Settings.hh"
#include "Achilles/System.hh"
#include "Achilles/Utilities.hh"
#include "fmt/ranges.h"
#include "spdlog/spdlog.h"
#include <functional>
#include <iostream>

// Demangling taken from: https://stackoverflow.com/a/34916852
#ifdef __GNUG__
#include <cxxabi.h>
std::string achilles::demangle(const char *name) {
    int status = -4;
    std::unique_ptr<char, void (*)(void *)> res{abi::__cxa_demangle(name, NULL, NULL, &status),
                                                std::free};

    return status == 0 ? res.get() : name;
}
#else
std::string achilles::demangle(const char *name) {
    return name;
}
#endif

using achilles::Settings;

Settings::Settings(const std::string &filename) {
    m_settings = IncludeFile(filename);
    CheckRequired();
}

void Settings::Print() const {
    std::cout << m_settings << std::endl;
}

const YAML::Node Settings::operator[](const std::string_view &key) const {
    auto paths = ParsePath(key);
    spdlog::trace("Settings: key = {}, parsed path = {}", key,
                  fmt::join(paths.begin(), paths.end(), ", "));
    if(paths.size() == 1) return m_settings[paths[0]];
    try {
        return seek(paths, m_settings);
    } catch(const SettingsError &) {
        auto msg = fmt::format("Settings: key {} can not be found", key);
        throw SettingsError(msg);
    }
}

YAML::Node Settings::operator[](const std::string_view &key) {
    auto paths = ParsePath(key);
    if(paths.size() == 1) return m_settings[paths[0]];
    spdlog::trace("Settings: key = {}, parsed path = {}", key,
                  fmt::join(paths.begin(), paths.end(), ", "));
    try {
        return seek(paths, m_settings);
    } catch(const SettingsError &) {
        auto msg = fmt::format("Settings: key {} can not be found", key);
        throw SettingsError(msg);
    }
}

bool Settings::Exists(const std::string_view &key) const {
    try {
        return static_cast<bool>((*this)[key]);
    } catch(const SettingsError &) { return false; }
}

YAML::Node Settings::IncludeFile(const std::string &filename) {
    YAML::Node node = YAML::LoadFile(Filesystem::FindFile(filename, "Settings"));

    MutableYAMLVisitor([](YAML::Node scalar) {
        if(scalar.Tag() == "!include") { scalar = IncludeFile(scalar.Scalar()); }
    })(node);

    return node;
}

void Settings::CheckRequired() const {
    for(const auto &option : m_required_options) {
        spdlog::trace("Looking for option {}", option);
        if(!(*this)[option]) {
            spdlog::debug("Printing out processed run card");
            if(spdlog::get_level() == spdlog::level::debug) Print();
            auto msg = fmt::format("Settings: Required option {} is not defined", option);
            throw SettingsError(msg);
        }
    }
}

std::vector<std::string> Settings::ParsePath(const std::string_view &keys) const {
    std::vector<std::string> tokens;
    tokenize(keys, tokens, "/");
    std::reverse(tokens.begin(), tokens.end());
    return tokens;
}

YAML::Node Settings::seek(std::vector<std::string> &keys, YAML::Node start) {
    if(!start.IsMap()) throw SettingsError("Settings: Invalid node to seek");
    auto key = keys.back();
    keys.pop_back();
    for(auto it = start.begin(); it != start.end(); ++it) {
        if(it->first.Scalar() == key) {
            if(keys.size() == 0) return it->second;
            return seek(keys, it->second);
        }
    }
    throw SettingsError("Key not found");
}
