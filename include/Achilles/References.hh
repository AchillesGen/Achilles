#ifndef REFERENCES_HH
#define REFERENCES_HH

#include <string>
#include <utility>
#include <vector>

#include "fmt/format.h"

namespace achilles {

using ref_map = std::vector<std::pair<std::string, std::string>>;

struct Reference {
    std::string m_type, m_label;
    ref_map m_data;

    Reference(std::string type, std::string label) 
        : m_type{std::move(type)}, m_label{std::move(label)} {}
    Reference(std::string type, std::string label, ref_map data)
        : m_type(std::move(type)), m_label{std::move(label)}, m_data{std::move(data)} {}

    void AddField(const std::string &key, const std::string &value) {
        m_data.emplace_back(key, value); 
    }

    std::string GetReference() const {
        std::string result;
        result += fmt::format("@{}{{},\n", m_type, m_label);
        for(const auto &entry : m_data) {
            result += fmt::format("\t{} = {},\n", entry.first, entry.second);
        }
        result += fmt::format("}");
        return result;
    }

    void WriteReference() const {
        fmt::print("@{}{{},\n", m_type, m_label);
        for(const auto &entry : m_data) {
            fmt::print("\t{} = {},\n", entry.first, entry.second);
        }
        fmt::print("}");
    };
};

}

#endif
