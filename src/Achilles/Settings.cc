#include "Achilles/Settings.hh"
#include "Achilles/System.hh"
#include "Achilles/Utilities.hh"
#include "fmt/ranges.h"
#include "spdlog/spdlog.h"
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

Rule::Rule(const std::string &_type, const YAML::Node &_options, const std::string &_description,
           const std::string &_error_message)
    : type(_type), description(_description), error_message(_error_message) {
    if(type.empty()) throw SettingsError("SettingsValidator: Rule type cannot be empty");
    if(description.empty())
        throw SettingsError("SettingsValidator: Rule description cannot be empty");
    if(error_message.empty())
        throw SettingsError("SettingsValidator: Rule error message cannot be empty");

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
            spdlog::debug("Loading rule {}: {}", rule["Type"].as<std::string>(),
                          rule["Description"].as<std::string>());
            auto type = rule["Type"].as<std::string>();
            auto options = rule;
            auto description = rule["Description"].as<std::string>();
            auto error_message = rule["ErrorMessage"].as<std::string>();

            AddRule(Rule(type, options, description, error_message));
        }
    } catch(const YAML::Exception &e) {
        spdlog::error("Settings: Failed to load rules from file {}: {}", filename, e.what());
        return false;
    }
    return true;
}

bool SettingsValidator::Validate(const Settings &settings) const {
    bool valid = true;
    for(const auto &rule : m_rules) {
        spdlog::debug("Validating rule: {}", rule.type);
        if(rule.type == "RequiredFields") {
            valid &= ValidateRequiredFields(settings, rule);
        } else if(rule.type == "ConditionalInteractionOption") {
            valid &= ValidateConditionalInteractionOption(settings, rule);
        } else if(rule.type == "RangeConstraint") {
            valid &= ValidateRangeConstraint(settings, rule);
        } else {
            spdlog::error("Settings: Unknown rule type: {}", rule.type);
            return false;
        }
    }
    return valid;
}

bool SettingsValidator::ValidateRequiredFields(const Settings &settings, const Rule &rule) const {
    for(const auto &option : rule.options.at("Fields")) {
        spdlog::trace("Checking required option: {}", option.as<std::string>());
        if(!settings.Exists(option.as<std::string>())) {
            spdlog::error(rule.error_message, option.as<std::string>());
            return false;
        }
    }
    return true;
}

bool SettingsValidator::ValidateConditionalInteractionOption(const Settings &settings,
                                                             const Rule &rule) const {
    if(settings.GetAs<bool>("Cascade/Run") == false) {
        spdlog::debug("Cascade is not running, skipping ConsistentInteractions validation");
        return true; // No interactions to validate if Cascade is not running
    }
    const auto &interactions = settings["Cascade/Interactions"];

    const auto &condition = rule.options.at("Condition");
    const auto &target = rule.options.at("Target");

    // Collect condition names
    std::vector<std::string> condition_names;
    if(condition["Name"]) {
        condition_names.push_back(condition["Name"].as<std::string>());
    } else if(condition["Interactions"]) {
        for(const auto &name : condition["Interactions"]) {
            condition_names.push_back(name.as<std::string>());
        }
    } else {
        spdlog::error("SettingsValidator: Condition must have a 'Name' or 'Interactions' field");
        return false;
    }
    const std::string target_name = target["Name"].as<std::string>();
    const YAML::Node &requirements = target["RequiredOptions"];

    bool condition_met = false;
    YAML::Node target_node;

    for(size_t i = 0; i < interactions.size(); ++i) {
        const YAML::Node &interaction = interactions[i];
        const std::string name = interaction["Name"].as<std::string>();

        if(std::find(condition_names.begin(), condition_names.end(), name) !=
           condition_names.end()) {
            condition_met = true;
            spdlog::debug("Found condition: {}", name);
        }
        if(name == target_name) {
            target_node = interaction["Options"];
            spdlog::debug("Found target: {}", name);
        }
    }

    if(!condition_met || !target_node) {
        return true; // If condition or target not found, no need to validate
    }

    for(const auto &vals : requirements) {
        const std::string key = vals.first.as<std::string>();
        const std::string value = vals.second.as<std::string>();

        if(!target_node[key]) {
            spdlog::error("Settings: Target {} does not have required option {}", target_name, key);
            return false;
        }

        const std::string target_value = target_node[key].as<std::string>();
        if(target_value != value) {
            spdlog::error("Settings: Target {} has invalid value for {}: expected {}, "
                          "got {}",
                          target_name, key, value, target_value);
            return false;
        }
    }
    return true;
}

bool SettingsValidator::ValidateRangeConstraint(const Settings &settings, const Rule &rule) const {
    if(rule.options.find("Path") != rule.options.end()) {
        return ValidateRangeConstraintScalar(settings, rule);
    } else if(rule.options.find("List") != rule.options.end()) {
        return ValidateRangeConstraintList(settings, rule);
    } else {
        spdlog::error("Settings: RangeConstraint rule must specify either 'Path' or 'List'");
        return false;
    }
}

bool SettingsValidator::ValidateRangeConstraintList(const Settings &settings,
                                                    const Rule &rule) const {
    const auto list_path = rule.options.at("List").as<std::string>();
    const auto param_path = rule.options.at("Parameter").as<std::string>();
    YAML::Node list_node;
    try {
        list_node = settings[list_path];
    } catch(const SettingsError &e) {
        if(rule.options.at("Optional").as<bool>()) {
            spdlog::debug("Settings: List {} is optional and does not exist", list_path);
            return true; // Optional list, no validation needed
        } else {
            throw e;
        }
    }
    if(!list_node || !list_node.IsSequence()) { return true; }
    spdlog::debug("Validating range constraint for list at path: {}", list_path);

    const auto &match = rule.options.at("Match");
    auto matches_conditions = [&](const YAML::Node &item,
                                  const YAML::Node &match_criteria) -> bool {
        for(const auto &kv : match_criteria) {
            std::string key = kv.first.as<std::string>();
            const YAML::Node &allowed = kv.second;
            if(!item[key]) return false;

            const std::string value = item[key].as<std::string>();

            if(allowed.IsSequence()) {
                bool found = false;
                for(const auto &option : allowed) {
                    if(value == option.as<std::string>()) {
                        found = true;
                        break;
                    }
                }
                if(!found) return false;
            } else if(allowed.IsScalar()) {
                // Check if value matches the single allowed scalar
                if(value != allowed.as<std::string>()) { return false; }
            } else {
                return false;
            }
        }
        return true;
    };

    for(size_t i = 0; i < list_node.size(); ++i) {
        const YAML::Node item = list_node[i];
        if(!matches_conditions(item, match)) continue;

        // Navigate to sub-parameter
        auto rel_path = settings.ParsePath(param_path);
        spdlog::debug("Validating parameter relative_path: {}", Settings::FormatPath(rel_path));
        auto value_nodes = Settings::seek(rel_path, item);

        if(value_nodes.empty()) {
            spdlog::error("Settings: Parameter {} not found in item at index {}", param_path, i);
            return false;
        }

        for(const auto &param_node : value_nodes) {
            if(!param_node.IsScalar()) {
                spdlog::error("Settings: Parameter {} at index {} is not a scalar", param_path, i);
                return false;
            }

            double value = param_node.as<double>();
            if(!ValidateRange(value, rule,
                              list_path + "[" + std::to_string(i) + "]/" + param_path)) {
                return false;
            }
        }
    }

    return true;
}

bool SettingsValidator::ValidateRangeConstraintScalar(const Settings &settings,
                                                      const Rule &rule) const {
    const auto &param_path = rule.options.at("Path").as<std::string>();
    spdlog::trace("Validating range constraint for path: {}", param_path);

    if(!settings.Exists(param_path)) {
        if(rule.options.at("Optional").as<bool>()) {
            spdlog::debug("Settings: Path {} is optional and does not exist", param_path);
            return true; // Optional path, no validation needed
        } else {
            spdlog::error("Settings: Path {} does not exist and is required", param_path);
            return false;
        }
    }

    double value = settings.GetAs<double>(param_path);
    return ValidateRange(value, rule, param_path);
}

bool SettingsValidator::ValidateRange(double value, const Rule &rule,
                                      const std::string &path) const {
    if(rule.options.find("Min") != rule.options.end()) {
        double min = rule.options.at("Min").as<double>();
        if(value < min) {
            spdlog::error("Settings: Value {} for path {} must be greater than or equal to {}",
                          value, path, min);
            return false;
        }
    }
    if(rule.options.find("Max") != rule.options.end()) {
        double max = rule.options.at("Max").as<double>();
        if(value > max) {
            spdlog::error("Settings: Value {} for path {} must be less than or equal to {}", value,
                          path, max);
            return false;
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
    if(!m_validator.Validate(*this)) { throw SettingsError("Settings: Validation failed"); }
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
    auto matches = seek(paths, m_settings);
    spdlog::trace("Settings: key = {}, parsed path = {}", key, FormatPath(paths));
    if(matches.empty()) {
        auto msg = fmt::format("Settings: key {} can not be found", key);
        throw SettingsError(msg);
    }
    if(matches.size() > 1) {
        spdlog::warn("Settings: Multiple results found for key {}, returning first result", key);
    }

    return matches.front();
}

YAML::Node Settings::operator[](const std::string_view &key) {
    auto paths = ParsePath(key);
    auto matches = seek(paths, m_settings);
    spdlog::trace("Settings: key = {}, parsed path = {}", key, FormatPath(paths));
    if(matches.empty()) {
        auto msg = fmt::format("Settings: key {} can not be found", key);
        throw SettingsError(msg);
    }
    if(matches.size() > 1) {
        spdlog::warn("Settings: Multiple results found for key {}, returning first result", key);
    }

    return matches.front();
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

std::vector<Settings::PathElement> Settings::ParsePath(const std::string_view &path) const {
    std::vector<Settings::PathElement> path_elements;
    std::stringstream ss((std::string(path)));
    std::string segment;

    while(std::getline(ss, segment, '/')) {
        if(segment.empty()) continue; // Skip empty segments
        if(segment == "*") {
            path_elements.emplace_back(WildcardTag{});
        } else if(std::all_of(segment.begin(), segment.end(), ::isdigit)) {
            path_elements.emplace_back(std::stoul(segment));
        } else {
            path_elements.push_back(segment);
        }
    }

    return path_elements;
}

std::vector<YAML::Node> Settings::seek(std::vector<PathElement> &path, YAML::Node start,
                                       size_t index) {
    if(!start) return {};
    if(index >= path.size()) return {start};

    const PathElement &element = path[index];
    std::vector<YAML::Node> result;

    if(std::holds_alternative<std::string>(element)) {
        const std::string &key = std::get<std::string>(element);
        YAML::Node next = start[key];
        if(next) {
            auto sub = seek(path, next, index + 1);
            result.insert(result.end(), sub.begin(), sub.end());
        }
    } else if(std::holds_alternative<size_t>(element)) {
        size_t idx = std::get<size_t>(element);
        if(!start.IsSequence() || idx >= start.size()) {
            throw SettingsError(
                fmt::format("Settings: Index {} out of bounds for path {}", idx, FormatPath(path)));
        }
        auto sub = seek(path, start[idx], index + 1);
        result.insert(result.end(), sub.begin(), sub.end());
    } else if(std::holds_alternative<WildcardTag>(element)) {
        if(start.IsSequence()) {
            for(size_t i = 0; i < start.size(); ++i) {
                auto sub = seek(path, start[i], index + 1);
                result.insert(result.end(), sub.begin(), sub.end());
            }
        } else if(start.IsMap()) {
            for(const auto &pair : start) {
                auto sub = seek(path, pair.second, index + 1);
                result.insert(result.end(), sub.begin(), sub.end());
            }
        } else {
            throw SettingsError(
                fmt::format("Settings: Wildcard path {} can only be used with sequences or maps",
                            FormatPath(path)));
        }
    }

    if(result.empty()) {
        throw SettingsError(fmt::format("Settings: Path {} not found", FormatPath(path)));
    }
    return result;
}
