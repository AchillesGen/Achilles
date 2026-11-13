// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef CONSOLE_HH
#define CONSOLE_HH

#include <cstddef>
#include <mutex>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace achilles::console {

struct Style {
    bool unicode{true};
    bool color{false};
    bool interactive{false};
};

// Detected on first use from ACHILLES_ASCII, TERM, NO_COLOR and whether stdout is a terminal
Style &GetStyle();
// Batch mode disables in-place redraws so logs stay readable
void Configure(bool batch_mode, bool force_ascii = false);

size_t DisplayWidth(std::string_view);
size_t TerminalWidth();
// 0 when not on a terminal
size_t RawWidth();

enum class Align { left, right, center };
std::string Pad(std::string_view, size_t width, Align align = Align::left);
std::string Truncate(std::string_view, size_t width);

std::string FormatDuration(double seconds);
std::string FormatCount(double count);
std::string FormatBytes(double bytes);
std::string FormatResult(double mean, double error);

std::vector<std::string> Wrap(const std::vector<std::string> &items, size_t width,
                              std::string_view separator = ", ");

enum class Color { none, bold, dim, red, green, yellow, blue, magenta, cyan };
// No-op unless color output is enabled
std::string Colorize(std::string_view, Color);

// "always" is for the final summary, which must survive quiet mode (-q)
enum class Verbosity { info, debug, always };

// Rows can be streamed with HeaderLines/RowLine/FooterLine for per-iteration output
class Table {
  public:
    enum class Style { grid, frame };

    explicit Table(std::vector<std::string> headers, std::vector<Align> aligns = {});

    void SetTitle(std::string title) { m_title = std::move(title); }
    void SetStyle(Style style) { m_style = style; }
    // 0 means the terminal width on a terminal and no limit otherwise
    void SetMaxWidth(size_t width) { m_max_width = width; }
    // Needed when streaming, since later rows are not known yet
    void SetMinWidths(std::vector<size_t> widths);
    void AddRow(std::vector<std::string> row);
    void AddSeparator();

    std::vector<std::string> Render() const;
    void Print(Verbosity verbosity = Verbosity::info) const;

    std::vector<std::string> HeaderLines() const;
    std::string RowLine(const std::vector<std::string> &row) const;
    std::string FooterLine() const;

  private:
    size_t MaxWidth() const;
    size_t InnerWidth(const std::vector<size_t> &) const;
    std::vector<size_t> Widths() const;
    std::string RuleLine(const char *left, const char *mid, const char *right,
                         const std::vector<size_t> &) const;
    std::string MiddleRule(const std::vector<size_t> &) const;
    std::vector<std::string> FormatRow(const std::vector<std::string> &,
                                       const std::vector<size_t> &, bool header) const;

    std::string m_title;
    Style m_style{Style::grid};
    size_t m_max_width{0};
    std::vector<std::string> m_headers;
    std::vector<Align> m_aligns;
    std::vector<size_t> m_min_widths;
    // An empty row marks a separator
    std::vector<std::vector<std::string>> m_rows;
};

// A row with an empty name continues the one above it
std::vector<std::string> Fields(const std::string &title,
                                const std::vector<std::pair<std::string, std::string>> &rows);

std::vector<std::string> Box(const std::vector<std::vector<std::string>> &sections,
                             size_t inner_width);

// No log prefix, but still reaches the log file once the loggers are set up
void Print(const std::vector<std::string> &lines, Verbosity verbosity = Verbosity::info);
void Print(const std::string &line, Verbosity verbosity = Verbosity::info);

// Log messages suspend it so they don't interleave with progress bar redraws
class LiveLine {
  public:
    static LiveLine &Instance();
    void Draw(const std::string &text);
    void Clear();
    void Suspend();
    void Resume();

  private:
    LiveLine() = default;
    void Write(const std::string &);
    std::recursive_mutex m_mutex;
    std::string m_text;
    bool m_drawn{false};
};

// Shares the main logger's sinks without the time and level prefix
inline constexpr const char *raw_logger_name = "achilles.raw";

void SetupLoggers(int level, int log_level, int flush_time, const std::string &destination);
// spdlog numbering
int ConsoleLevel();

} // namespace achilles::console

#endif
