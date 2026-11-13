// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Contributors.hh"
#include "Achilles/ContributorsData.hh"

#include "fmt/format.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <stdexcept>

using achilles::Contributors;

namespace {

std::string ParseName(const YAML::Node &node, const std::string &section) {
    if(!node.IsScalar())
        throw std::runtime_error(
            fmt::format("Contributors: entries in \"{}\" must be plain names", section));
    auto name = node.as<std::string>();
    if(name.empty())
        throw std::runtime_error(fmt::format("Contributors: empty name in \"{}\"", section));
    return name;
}

// Sort key: the surname (last word), then the full name, ignoring case
std::string SortKey(const std::string &name) {
    const auto last = name.find_last_of(' ');
    std::string key = (last == std::string::npos ? name : name.substr(last + 1)) + " " + name;
    std::transform(key.begin(), key.end(), key.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return key;
}

std::vector<std::string> ParseSection(const YAML::Node &root, const std::string &section) {
    std::vector<std::string> people;
    const auto node = root[section];
    if(!node || node.IsNull()) return people;
    if(!node.IsSequence())
        throw std::runtime_error(
            fmt::format("Contributors: \"{}\" must be a list of names", section));
    for(const auto &entry : node) people.push_back(ParseName(entry, section));
    std::stable_sort(people.begin(), people.end(),
                     [](const auto &a, const auto &b) { return SortKey(a) < SortKey(b); });
    return people;
}

} // namespace

Contributors Contributors::Parse(const std::string &yaml) {
    YAML::Node root;
    try {
        root = YAML::Load(yaml);
    } catch(const YAML::Exception &error) {
        throw std::runtime_error(fmt::format("Contributors: invalid YAML: {}", error.what()));
    }
    if(!root.IsMap()) throw std::runtime_error("Contributors: expected a map of sections");

    for(const auto &section : root) {
        const auto key = section.first.as<std::string>();
        if(key != "authors" && key != "former_authors" && key != "undergraduate_contributors")
            throw std::runtime_error(fmt::format("Contributors: unknown section \"{}\"", key));
    }

    Contributors result;
    result.authors = ParseSection(root, "authors");
    result.former_authors = ParseSection(root, "former_authors");
    result.undergraduate_contributors = ParseSection(root, "undergraduate_contributors");
    return result;
}

const std::string &Contributors::BuiltinYAML() {
    static const std::string yaml{generated::contributors_yaml};
    return yaml;
}

const Contributors &Contributors::Builtin() {
    static const Contributors contributors = [] {
        try {
            return Parse(BuiltinYAML());
        } catch(const std::runtime_error &error) {
            // The banner is not worth failing a run over; the test suite checks the file
            std::fprintf(stderr, "%s\n", error.what());
            return Contributors{};
        }
    }();
    return contributors;
}
