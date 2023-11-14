#pragma once

#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include "fmt/core.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

std::string demangle(const char *name);

class SettingsError : public std::runtime_error {
    public:
        SettingsError(const std::string &what) : std::runtime_error(what) {}
};

template<bool is_const>
class YAMLVisitor {
    public:
        using value_type = YAML::Node;
        using cv_type = std::conditional_t<is_const, const YAML::Node, YAML::Node>;
        using iterator_type = std::conditional_t<is_const, YAML::const_iterator, YAML::iterator>;
        using Callback = std::function<void(cv_type)>;
        YAMLVisitor(Callback callback) : m_seen{}, m_callback{callback} {}

        void operator()(cv_type root) {
            m_seen.push_back(root);
            if(root.IsMap()) {
                for(iterator_type it = root.begin(); it != root.end(); ++it) {
                    next(it->second);
                }
            } else if(root.IsSequence()) {
                for(iterator_type it = root.begin(); it != root.end(); ++it) {
                    next(*it);
                }
            } else if(root.IsScalar()) {
                m_callback(root);
            }
        }

    private:
        void next(cv_type target) {
            if(std::find(m_seen.begin(), m_seen.end(), target) == m_seen.end())
                (*this)(target);
        }

        std::vector<value_type> m_seen;
        Callback m_callback;
};

using MutableYAMLVisitor = YAMLVisitor<false>;
using ConstYAMLVisitor = YAMLVisitor<true>;

class Settings {
    public:
        Settings(const std::string &filename);
        static Settings& MainSettings() {
            return main_settings;
        }
        static void LoadMainSettings(const std::string &filename) {
            main_settings = Settings(filename);
        }
        void Print() const;

        template<typename type>
        type GetAs(const std::string_view &key) const {
            try {
                return (*this)[key].as<type>();
            } catch (const YAML::TypedBadConversion<type> &e) {
                auto msg = fmt::format("Settings: Option {} is not of type {}",
                                       key, demangle(typeid(type).name()));
                throw SettingsError(msg);
            }
        }

        YAML::Node Root() const { return m_settings; }
        const YAML::Node operator[](const std::string_view &key) const;
        YAML::Node operator[](const std::string_view &key);
        bool Exists(const std::string_view &key) const;

    private:
        static YAML::Node IncludeFile(const std::string &filename);
        void CheckRequired() const;
        std::vector<std::string> ParsePath(const std::string_view &) const;
        static YAML::Node seek(std::vector<std::string> &keys, YAML::Node start);

        YAML::Node m_settings;
        static constexpr std::array<std::string_view, 9> m_required_options = {
            "Main/NEvents", "Main/Output/Format", "Main/Output/Name", "Process/Final States",
            "NuclearModel/Model", "Nucleus", "Cascade/Run", "Options/Unweighting", "Options/Initialize" };
        static Settings main_settings;
};

}
