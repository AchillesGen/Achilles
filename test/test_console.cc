// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Console.hh"
#include "Achilles/Contributors.hh"
#include "Achilles/Splash.hh"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_string.hpp"

#include <stdexcept>

using namespace achilles::console;

namespace {

// Sets the output style for the duration of a test
class ScopedStyle {
  public:
    ScopedStyle(bool unicode, bool color) : m_saved{GetStyle()} {
        GetStyle().unicode = unicode;
        GetStyle().color = color;
        GetStyle().interactive = false;
    }
    ~ScopedStyle() { GetStyle() = m_saved; }

  private:
    Style m_saved;
};

} // namespace

TEST_CASE("Display width ignores escapes and counts code points", "[console]") {
    CHECK(DisplayWidth("abc") == 3);
    CHECK(DisplayWidth("Gutiérrez") == 9);
    CHECK(DisplayWidth("╭─╮") == 3);
    CHECK(DisplayWidth("\033[1mbold\033[0m") == 4);
}

TEST_CASE("Padding and truncation", "[console]") {
    ScopedStyle style(true, false);
    CHECK(Pad("ab", 4) == "ab  ");
    CHECK(Pad("ab", 4, Align::right) == "  ab");
    CHECK(Pad("ab", 5, Align::center) == " ab  ");
    CHECK(Pad("toolong", 3) == "toolong");
    CHECK(Pad("é", 3) == "é  ");

    CHECK(Truncate("short", 10) == "short");
    CHECK(Truncate("abcdefgh", 5) == "abcd…");
    CHECK(DisplayWidth(Truncate("\033[1mabcdefgh\033[0m", 5)) == 5);
}

TEST_CASE("Durations and counts are human readable", "[console]") {
    CHECK(FormatDuration(0.25) == "250ms");
    CHECK(FormatDuration(4.25) == "4.2s");
    CHECK(FormatDuration(42) == "42s");
    CHECK(FormatDuration(245) == "4m 05s");
    CHECK(FormatDuration(2 * 3600 + 3 * 60 + 10) == "2h 03m");
    CHECK(FormatDuration(86400 + 4 * 3600) == "1d 04h");
    CHECK(FormatDuration(-1) == "--");

    CHECK(FormatCount(950) == "950");
    CHECK(FormatCount(12345) == "12.3k");
    CHECK(FormatCount(100000) == "100k");
    CHECK(FormatCount(4560000) == "4.56M");
}

TEST_CASE("Wrap joins items without exceeding the width", "[console]") {
    const std::vector<std::string> names{"Alpha Beta", "Gamma Delta", "Epsilon", "Zeta Eta"};
    const auto lines = Wrap(names, 25);
    REQUIRE(lines.size() == 2);
    CHECK(lines[0] == "Alpha Beta, Gamma Delta,");
    CHECK(lines[1] == "Epsilon, Zeta Eta");
    for(const auto &line : lines) CHECK(DisplayWidth(line) <= 25);

    CHECK(Wrap({}, 10).empty());
    CHECK(Wrap({"Single"}, 3) == std::vector<std::string>{"Single"});
}

TEST_CASE("Tables are drawn with box characters", "[console]") {
    Table table({"Name", "Value"}, {Align::left, Align::right});
    table.AddRow({"a", "1"});
    table.AddSeparator();
    table.AddRow({"long name", "22"});

    SECTION("Unicode") {
        ScopedStyle style(true, false);
        const std::vector<std::string> expected{
            "╭───────────┬───────╮", "│ Name      │ Value │", "├───────────┼───────┤",
            "│ a         │     1 │", "├───────────┼───────┤", "│ long name │    22 │",
            "╰───────────┴───────╯",
        };
        CHECK(table.Render() == expected);
    }

    SECTION("ASCII") {
        ScopedStyle style(false, false);
        const auto lines = table.Render();
        CHECK(lines.front() == "+-----------+-------+");
        CHECK(lines[1] == "| Name      | Value |");
    }

    SECTION("Title widens the table when needed") {
        ScopedStyle style(true, false);
        table.SetTitle("A rather long table title");
        const auto lines = table.Render();
        for(const auto &line : lines) CHECK(DisplayWidth(line) == DisplayWidth(lines[0]));
        CHECK_THAT(lines[1], Catch::Matchers::ContainsSubstring("A rather long table title"));
        CHECK(lines[2] == "├───────────┬───────────────┤");
    }

    SECTION("Frame style draws only the outer border") {
        ScopedStyle style(true, false);
        table.SetStyle(Table::Style::frame);
        const std::vector<std::string> expected{
            "╭──────────────────╮", "│ Name       Value │", "├──────────────────┤",
            "│ a              1 │", "├──────────────────┤", "│ long name     22 │",
            "╰──────────────────╯",
        };
        CHECK(table.Render() == expected);
    }

    SECTION("Tables wider than the limit wrap their widest column") {
        ScopedStyle style(true, false);
        table.AddRow({"a much longer name than fits", "3"});
        table.SetMaxWidth(20);
        const auto lines = table.Render();
        for(const auto &line : lines) CHECK(DisplayWidth(line) == 20);
        CHECK(lines.size() > 9);
    }

    SECTION("Streamed rows line up with the header") {
        ScopedStyle style(true, false);
        Table streamed({"Iter", "Result"}, {Align::right, Align::right});
        streamed.SetMinWidths({4, 8});
        const auto header = streamed.HeaderLines();
        const auto row = streamed.RowLine({"1", "3.14"});
        CHECK(DisplayWidth(row) == DisplayWidth(header.front()));
        CHECK(DisplayWidth(streamed.FooterLine()) == DisplayWidth(header.front()));
    }
}

TEST_CASE("Boxes pad every line to the same width", "[console]") {
    ScopedStyle style(true, false);
    const auto lines = Box({{"one", "two"}, {"three"}}, 10);
    REQUIRE(lines.size() == 6);
    CHECK(lines[0] == "╔══════════╗");
    CHECK(lines[1] == "║one       ║");
    CHECK(lines[3] == "╟──────────╢");
    CHECK(lines[5] == "╚══════════╝");
}

TEST_CASE("Contributors are read from YAML", "[contributors]") {
    using achilles::Contributors;

    SECTION("Sections are lists of names") {
        const auto people = Contributors::Parse(R"(
authors:
  - Jane Doe
  - John Roe
former_authors: []
undergraduate_contributors:
  - Sam Poe
)");
        REQUIRE(people.authors.size() == 2);
        CHECK(people.authors[0] == "Jane Doe");
        CHECK(people.authors[1] == "John Roe");
        CHECK(people.former_authors.empty());
        REQUIRE(people.undergraduate_contributors.size() == 1);
        CHECK(people.undergraduate_contributors[0] == "Sam Poe");
    }

    SECTION("Names are sorted by surname") {
        const auto people = Contributors::Parse("authors: [Jane Roe, john doe, Ann Doe, Zed]");
        REQUIRE(people.authors.size() == 4);
        CHECK(people.authors[0] == "Ann Doe");
        CHECK(people.authors[1] == "john doe");
        CHECK(people.authors[2] == "Jane Roe");
        CHECK(people.authors[3] == "Zed");
    }

    SECTION("Missing sections are empty") {
        const auto people = Contributors::Parse("authors: [Jane Doe]");
        CHECK(people.authors.size() == 1);
        CHECK(people.former_authors.empty());
        CHECK(people.undergraduate_contributors.empty());
    }

    SECTION("Malformed input is rejected") {
        CHECK_THROWS_AS(Contributors::Parse("authors: Jane Doe"), std::runtime_error);
        CHECK_THROWS_AS(Contributors::Parse("authors: [{name: Jane Doe}]"), std::runtime_error);
        CHECK_THROWS_AS(Contributors::Parse("authors: [[Jane Doe]]"), std::runtime_error);
        CHECK_THROWS_AS(Contributors::Parse("autors: [Jane Doe]"), std::runtime_error);
        CHECK_THROWS_AS(Contributors::Parse("authors: [Jane"), std::runtime_error);
    }

    SECTION("The embedded CONTRIBUTORS.yml is valid") {
        const auto people = Contributors::Parse(Contributors::BuiltinYAML());
        CHECK_FALSE(people.authors.empty());
    }
}

TEST_CASE("Splash lines all have the same width", "[contributors]") {
    ScopedStyle style(true, false);
    const auto lines = achilles::SplashLines("1.2.3");
    REQUIRE_FALSE(lines.empty());
    for(const auto &line : lines) CHECK(DisplayWidth(line) == DisplayWidth(lines.front()));

    bool has_version = false;
    for(const auto &line : lines) has_version |= line.find("1.2.3") != std::string::npos;
    CHECK(has_version);
}
