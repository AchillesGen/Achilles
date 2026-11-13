// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Splash.hh"
#include "Achilles/Console.hh"
#include "Achilles/Contributors.hh"

#include <algorithm>
#include <string_view>

using achilles::console::Color;
using achilles::console::Colorize;

namespace {

constexpr size_t inner_width = 69;
constexpr size_t heading_indent = 4;
constexpr size_t name_indent = 8;

// clang-format off
const std::vector<std::string> title_art = {
    R"art(    .d8b.   .o88b. db   db d888888b db      db      d88888b .d8888.)art",
    R"art(   d8' `8b d8P  Y8 88   88   `88'   88      88      88'     88'  YP)art",
    R"art(   88ooo88 8P      88ooo88    88    88      88      88ooooo `8bo.)art",
    R"art(   88~~~88 8b      88~~~88    88    88      88      88~~~~~  `Y8b.)art",
    R"art(   88   88 Y8b  d8 88   88   .88.   88booo. 88booo. 88.     db  8D)art",
    R"art(   YP   YP  `Y88P' YP   YP Y888888P Y88888P Y88888P Y88888P `8888Y')art",
};

// Traced from assets/logo.svg: the helmet over the four electron orbits
const std::vector<std::string> emblem_art = {
    R"art(                            /----\)art",
    R"art(                          // .db. \\)art",
    R"art(                         / .d8888b. \)art",
    R"art(                        /  88888888' \)art",
    R"art(          ______       |    Y88888'   |       ______)art",
    R"art(       //       ----\\|      Y888'     |//----      \\)art",
    R"art(      |               |\\\   d88b   ///|              |)art",
    R"art(      |              |  .oo88888888bo.  |             |)art",
    R"art(       |             |o8888888888888888o|             |)art",
    R"art(       |\            d888888888888888888b.           /|)art",
    R"art(        \\  ____----d888888888888888888888----____  //)art",
    R"art(      _/\\-      .d8888888888888888888888b.      -//\_)art",
    R"art(   ////    \      d888888888888888888888888b      /    \\\\)art",
    R"art( //         \\    d888888888888888888888888b    //         \\)art",
    R"art(|/            \ //`88b ""Y8888888888P"" d88b\\//            \|)art",
    R"art(|             /\\  88b       `88b       d88' /\\             |)art",
    R"art( \\          /   \\888bo.    d88b    .od888b/   \          //)art",
    R"art(   \\\     //     d8888888b  d88b  d8888888b     \\     ///)art",
    R"art(      \\__/       d88888888. 8888  88888888b       \__//)art",
    R"art(         / ----___"88888888b "88"  88888888P___---- \)art",
    R"art(        /          `Y888888b\----/d888888P'          \)art",
    R"art(       |             `88888b \\// d88888"             |)art",
    R"art(       |             | "888b//  \\d888P |             |)art",
    R"art(      |               |/`Y8b      d8P'\|              |)art",
    R"art(       \           _//|   `b      d"   |\\_           /)art",
    R"art(        \\---------    |              |    ---------//)art",
    R"art(                        \            /)art",
    R"art(                         \          /)art",
    R"art(                          \        /)art",
    R"art(                           \\____//)art",
};
// clang-format on

void AddNames(std::vector<std::string> &lines, const std::string &heading,
              const std::vector<std::string> &names) {
    if(names.empty()) return;

    lines.push_back(std::string(heading_indent, ' ') + Colorize(heading, Color::bold));
    for(const auto &line : achilles::console::Wrap(names, inner_width - name_indent - 2))
        lines.push_back(std::string(name_indent, ' ') + line);
    lines.emplace_back();
}

// Orbits in blue and the helmet in bold, like the logo
std::string ColorEmblem(const std::string &line) {
    std::string result, run;
    bool orbit = false;
    auto flush = [&]() {
        if(run.empty()) return;
        result += Colorize(run, orbit ? Color::blue : Color::bold);
        run.clear();
    };
    for(const char c : line) {
        if(c == ' ') {
            flush();
            result += c;
            continue;
        }
        const bool is_orbit = std::string_view("-_/\\|").find(c) != std::string_view::npos;
        if(is_orbit != orbit) flush();
        orbit = is_orbit;
        run += c;
    }
    flush();
    return result;
}

// Indent that centers a block of lines as a whole, so that art keeps its shape
std::string BlockIndent(const std::vector<std::string> &lines) {
    size_t width = 0, margin = std::string::npos;
    for(const auto &line : lines) {
        width = std::max(width, line.size());
        const size_t lead = line.find_first_not_of(' ');
        if(lead != std::string::npos) margin = std::min(margin, lead);
    }
    if(margin == std::string::npos) margin = 0;
    // Leading spaces shared by every line are part of the margin, not the art
    const size_t art_width = width - margin;
    const size_t indent = art_width < inner_width ? (inner_width - art_width) / 2 : 0;
    return std::string(indent > margin ? indent - margin : 0, ' ');
}

} // namespace

std::vector<std::string> achilles::SplashLines(const std::string &version) {
    std::vector<std::string> title{""};
    const std::string title_indent = BlockIndent(title_art);
    for(const auto &line : title_art) title.push_back(title_indent + Colorize(line, Color::cyan));
    title.emplace_back();

    std::vector<std::string> emblem{""};
    const std::string emblem_indent = BlockIndent(emblem_art);
    for(const auto &line : emblem_art) emblem.push_back(emblem_indent + ColorEmblem(line));
    emblem.emplace_back();

    const auto &people = Contributors::Builtin();
    std::vector<std::string> credits{""};
    credits.push_back(std::string(heading_indent, ' ') + Colorize("Version", Color::bold) + "  " +
                      version);
    credits.emplace_back();
    AddNames(credits, "Authors", people.authors);
    AddNames(credits, "Former Authors", people.former_authors);
    AddNames(credits, "Undergraduate Contributors", people.undergraduate_contributors);

    return console::Box({title, emblem, credits}, inner_width);
}

void achilles::PrintSplash(const std::string &version) {
    console::Print(SplashLines(version));
}
