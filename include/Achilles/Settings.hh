#pragma once

#include "fmt/core.h"
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

std::string demangle(const char *name);

class SettingsError : public std::runtime_error {
  public:
    SettingsError(const std::string &what) : std::runtime_error(what) {}
};

template <bool is_const> class YAMLVisitor {
  public:
    using value_type = YAML::Node;
    using cv_type = std::conditional_t<is_const, const YAML::Node, YAML::Node>;
    using iterator_type = std::conditional_t<is_const, YAML::const_iterator, YAML::iterator>;
    using Callback = std::function<void(cv_type)>;
    YAMLVisitor(Callback callback) : m_seen{}, m_callback{callback} {}

    void operator()(cv_type root) {
        m_seen.push_back(root);
        if(root.IsMap()) {
            for(iterator_type it = root.begin(); it != root.end(); ++it) { next(it->second); }
        } else if(root.IsSequence()) {
            for(iterator_type it = root.begin(); it != root.end(); ++it) { next(*it); }
        } else if(root.IsScalar()) {
            m_callback(root);
        }
    }

  private:
    void next(cv_type target) {
        if(std::find(m_seen.begin(), m_seen.end(), target) == m_seen.end()) (*this)(target);
    }

    std::vector<value_type> m_seen;
    Callback m_callback;
};

using MutableYAMLVisitor = YAMLVisitor<false>;
using ConstYAMLVisitor = YAMLVisitor<true>;

class Settings;

class SettingsValidator {
  public:
    struct Rule {
        std::string name;
        std::string type;
        std::string description;
        std::string error_message;
        std::map<std::string, YAML::Node> options;

        Rule() = default;
        Rule(const std::string &_type, const YAML::Node &_options,
             const std::string &_description, const std::string &_error_message);
    };

    bool LoadRules(const std::string &filename);
    bool Validate(const Settings &settings) const;
    void AddRule(const Rule &rule) { m_rules.push_back(rule); }

  private:
    bool ValidateRequiredFields(const Settings &settings, const Rule &rule) const;
    bool ValidateConditionalInteractionOption(const Settings &settings, const Rule &rule) const;
    bool ValidateRangeConstraint(const Settings &settings, const Rule &rule) const;
    bool ValidateRangeConstraintList(const Settings &settings, const Rule &rule) const;
    bool ValidateRangeConstraintScalar(const Settings &settings, const Rule &rule) const;
    bool ValidateRange(double, const Rule &rule, const std::string &path) const;

    std::vector<Rule> m_rules;
};

class Settings {
  public:
    Settings(const std::string &filename, const std::string &rules="data/Rules.yml");
    Settings(const YAML::Node &node);
    static Settings &MainSettings() { return main_settings; }
    static void LoadMainSettings(const std::string &filename) {
        main_settings = Settings(filename);
    }
    void Print() const;

    template <typename type> type GetAs(const std::string_view &key) const {
        try {
            return (*this)[key].as<type>();
        } catch(const YAML::TypedBadConversion<type> &e) {
            auto msg = fmt::format("Settings: Option {} is not of type {}", key,
                                   demangle(typeid(type).name()));
            throw SettingsError(msg);
        }
    }

    YAML::Node Root() const { return m_settings; }
    const YAML::Node operator[](const std::string_view &key) const;
    YAML::Node operator[](const std::string_view &key);
    bool Exists(const std::string_view &key) const;

    struct WildcardTag {};
    using PathElement = std::variant<std::string, size_t, WildcardTag>;

    static std::string FormatPath(const std::vector<PathElement> &path) {
        std::string formatted_path;
        for(const auto &element : path) {
            if(!formatted_path.empty()) {
                formatted_path += "/";
            }
            if(std::holds_alternative<std::string>(element)) {
                formatted_path += std::get<std::string>(element);
            } else if(std::holds_alternative<size_t>(element)) {
                formatted_path += std::to_string(std::get<size_t>(element));
            } else {
                formatted_path += "*";
            }
        }
        return formatted_path;
    };

  private:
    static YAML::Node IncludeFile(const std::string &filename);
    std::vector<PathElement> ParsePath(const std::string_view &) const;
    static std::vector<YAML::Node> seek(std::vector<PathElement> &keys, YAML::Node start, size_t index=0);

    YAML::Node m_settings;
    SettingsValidator m_validator;
    static Settings main_settings;
    friend class SettingsValidator;
};

} // namespace achilles
