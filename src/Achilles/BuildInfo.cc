// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/BuildInfo.hh"
#include "Achilles/BuildConfig.hh"
#include "Achilles/Console.hh"
#include "Achilles/System.hh"
#include "Achilles/Version.hh"

#include "fmt/format.h"
#include "fmt/ranges.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <string>
#include <string_view>
#include <sys/utsname.h>

using achilles::InfoRow;

namespace {

bool Enabled(std::string_view value) {
    std::string upper(value);
    std::transform(upper.begin(), upper.end(), upper.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return upper == "ON" || upper == "TRUE" || upper == "1" || upper == "YES";
}

std::string CxxCompiler() {
#if defined(__clang__)
    return fmt::format("Clang {}", __clang_version__);
#elif defined(__GNUC__)
    return fmt::format("GCC {}", __VERSION__);
#else
    return "unknown";
#endif
}

std::string CxxStandard() {
    switch(__cplusplus) {
    case 201703L:
        return "C++17";
    case 202002L:
        return "C++20";
    case 202302L:
        return "C++23";
    default:
        return fmt::format("C++ ({})", __cplusplus);
    }
}

} // namespace

std::vector<InfoRow> achilles::BuildRows() {
    namespace gen = achilles::generated;
    std::vector<std::string> on, off;
    const std::vector<std::pair<std::string, std::string_view>> features{
        {"gzip", gen::feature_gzip},
        {"HepMC3/NuHepMC", gen::feature_hepmc3},
        {"Sherpa", gen::feature_sherpa},
        {"ROOT", gen::feature_root},
        {"Python", gen::feature_python},
        {"autodiff", gen::feature_autodiff},
        {"low memory", gen::feature_low_memory},
    };
    for(const auto &[name, value] : features) (Enabled(value) ? on : off).push_back(name);

    const std::string build_type =
        gen::build_type.empty() ? "not set" : std::string(gen::build_type);
    return {
        {"Version", ACHILLES_VERSION},
        {"Build type", build_type},
        {"C++ compiler", fmt::format("{} ({})", CxxCompiler(), CxxStandard())},
        {"Fortran compiler", std::string(gen::fortran_compiler)},
        {"CMake", std::string(gen::cmake_version)},
        {"Features on", on.empty() ? "none" : fmt::format("{}", fmt::join(on, ", "))},
        {"Features off", off.empty() ? "none" : fmt::format("{}", fmt::join(off, ", "))},
        {"Data release", std::string(gen::data_release)},
    };
}

std::vector<InfoRow> achilles::SystemRows() {
    std::vector<InfoRow> rows;
    utsname system{};
    if(uname(&system) == 0)
        rows.push_back(
            {"System", fmt::format("{} {} ({})", system.sysname, system.release, system.machine)});

    std::vector<std::string> dirs;
    for(const auto &dir : Filesystem::AchillesPath()) dirs.push_back(dir.string());
    rows.push_back({"Data search path", fmt::format("{}", fmt::join(dirs, ":"))});
    for(const char *name : {"ACHILLES_PATH", "ACHILLES_DATA_DIR", "ACHILLES_PLUGIN_PATH"}) {
        if(const char *value = std::getenv(name)) rows.push_back({name, value});
    }
    return rows;
}

std::string achilles::InfoText(const std::vector<InfoRow> &rows) {
    size_t width = 0;
    for(const auto &row : rows) width = std::max(width, row.name.size());
    std::string text;
    // A row without a name continues the one above it
    for(const auto &row : rows)
        text += fmt::format("  {:<{}}  {}\n", row.name.empty() ? "" : row.name + ":", width + 1,
                            row.value);
    return text;
}

void achilles::PrintInfo(const std::string &title, const std::vector<InfoRow> &rows, bool debug) {
    std::vector<std::pair<std::string, std::string>> fields;
    for(const auto &row : rows) fields.emplace_back(row.name, row.value);
    console::Print(console::Fields(title, fields),
                   debug ? console::Verbosity::debug : console::Verbosity::info);
}
