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

using achilles::SettingsValidator;
using Rule = achilles::SettingsValidator::Rule;
using achilles::Settings;

Rule::Rule(const std::string &_type, const YAML::Node &_options,
           const std::string &_description, const std::string &_error_message)
           : type(_type), description(_description), error_message(_error_message) {
    if(type.empty())
        throw SettingsError("SettingsValidator: Rule type cannot be empty");
    if(description.empty())
        throw SettingsError("SettingsValidator: Rule description cannot be empty");
    if(error_message.empty())
        throw SettingsError("SettingsValidator: Rule error message cannot be empty");

    std::cout << "options: " << _options << std::endl;
    for(const auto &option : _options) {
        auto key = option.first.as<std::string>();
        if(key == "Type" || key == "Description" || key == "ErrorMessage") {
            continue; // Skip reserved keys
        }
        options[key] = option.second;
    }
}


bool SettingsValidator::LoadRules(const std::string &filename) {
    try {
        YAML::Node rules = YAML::LoadFile(Filesystem::FindFile(filename, "Settings"));
        for(const auto &rule : rules["Rules"]) {
            std::cout << rule << std::endl;
            spdlog::debug("Loading rule {}: {}",
                          rule["Type"].as<std::string>(),
                          rule["Description"].as<std::string>());
            auto type = rule["Type"].as<std::string>();
            auto options = rule;
            auto description = rule["Description"].as<std::string>();
            auto error_message = rule["ErrorMessage"].as<std::string>();

            AddRule(Rule(type, options, description, error_message));
        }
    } catch(const YAML::Exception &e) {
        spdlog::error("SettingsValidator: Failed to load rules from file {}: {}", filename,
                      e.what());
        return false;
    }
    return true;
}

bool SettingsValidator::Validate(const Settings &settings) const {
    for(const auto &rule : m_rules) {
        spdlog::debug("Validating rule: {}", rule.type);
        if(rule.type == "RequiredFields") {
            for(const auto &option : rule.options.at("Fields")) {
                spdlog::trace("Checking global option: {}", option.as<std::string>());
                if(!settings.Exists(option.as<std::string>())) {
                    spdlog::error(rule.error_message, option.as<std::string>());
                    return false;
                }
            }
        }
    }
    return true;
}

Settings::Settings(const std::string &filename, const std::string &rules) {
    spdlog::trace("Loading settings from file: {}", filename);
    m_settings = IncludeFile(filename);
    if(!m_validator.LoadRules(rules)) {
        throw SettingsError("Settings: Failed to load validation rules from " + rules);
    }
    spdlog::trace("Settings loaded, validating...");
    if(!m_validator.Validate(*this)) {
        throw SettingsError("Settings: Validation failed");
    }
    // CheckRequired();
}

Settings::Settings(const YAML::Node &node) {
    m_settings = node;
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
