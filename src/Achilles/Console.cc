// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Console.hh"

#include "fmt/format.h"
#include "spdlog/pattern_formatter.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <sys/ioctl.h>
#include <unistd.h>

using namespace achilles::console;

namespace {

bool EnvSet(const char *name) {
    const char *value = std::getenv(name);
    return value != nullptr && value[0] != '\0';
}

Style DetectStyle() {
    Style style;
    const bool tty = isatty(fileno(stdout)) != 0;
    const char *term = std::getenv("TERM");
    const bool dumb = term != nullptr && std::string_view(term) == "dumb";

    style.unicode = !EnvSet("ACHILLES_ASCII") && !dumb;
    style.color = tty && !dumb && !EnvSet("NO_COLOR");
    style.interactive = tty && !dumb;
    return style;
}

struct Glyphs {
    const char *horizontal, *vertical;
    const char *top_left, *top_mid, *top_right;
    const char *mid_left, *mid_mid, *mid_right;
    const char *bottom_left, *bottom_mid, *bottom_right;
};

const Glyphs &TableGlyphs() {
    static const Glyphs unicode{"─", "│", "╭", "┬", "╮", "├", "┼", "┤", "╰", "┴", "╯"};
    static const Glyphs ascii{"-", "|", "+", "+", "+", "+", "+", "+", "+", "+", "+"};
    return GetStyle().unicode ? unicode : ascii;
}

const Glyphs &BoxGlyphs() {
    // The mid glyphs are used for the thin rule between sections
    static const Glyphs unicode{"═", "║", "╔", "═", "╗", "╟", "─", "╢", "╚", "═", "╝"};
    static const Glyphs ascii{"=", "|", "+", "=", "+", "+", "-", "+", "+", "=", "+"};
    return GetStyle().unicode ? unicode : ascii;
}

std::string Repeat(const char *glyph, size_t count) {
    std::string result;
    for(size_t i = 0; i < count; ++i) result += glyph;
    return result;
}

// Raw lines are indented to line up with the text of log messages ("12:00:00 INFO  ")
constexpr size_t level_width = 6, time_width = 9;

bool IsContinuationByte(char c) {
    return (static_cast<unsigned char>(c) & 0xC0) == 0x80;
}

std::string StripAnsi(std::string_view text) {
    std::string result;
    result.reserve(text.size());
    for(size_t i = 0; i < text.size(); ++i) {
        if(text[i] == '\033' && i + 1 < text.size() && text[i + 1] == '[') {
            i += 2;
            while(i < text.size() && (text[i] < '@' || text[i] > '~')) ++i;
            continue;
        }
        result += text[i];
    }
    return result;
}

} // namespace

size_t achilles::console::RawWidth() {
    if(!GetStyle().interactive) return 0;
    // Raw lines are indented, and a line that fills the last column wraps on some terminals
    const size_t indent = level_width + 1;
    const size_t width = TerminalWidth();
    return width > indent ? width - indent : 1;
}

namespace {

// A word longer than width is cut and loses its colors
std::vector<std::string> WrapText(const std::string &text, size_t width) {
    if(width == 0 || DisplayWidth(text) <= width) return {text};
    std::vector<std::string> lines;
    std::string line, word;
    auto add_word = [&]() {
        while(DisplayWidth(word) > width) {
            if(!line.empty()) lines.push_back(std::move(line));
            line.clear();
            const std::string plain = StripAnsi(word);
            size_t used = 0, cut = 0;
            for(; cut < plain.size(); ++cut) {
                if(!IsContinuationByte(plain[cut]) && used++ == width) break;
            }
            lines.push_back(plain.substr(0, cut));
            word = plain.substr(cut);
        }
        if(word.empty()) return;
        if(line.empty())
            line = word;
        else if(DisplayWidth(line) + 1 + DisplayWidth(word) <= width)
            line += " " + word;
        else {
            lines.push_back(std::move(line));
            line = word;
        }
        word.clear();
    };
    for(const char c : text) {
        if(c == ' ')
            add_word();
        else
            word += c;
    }
    add_word();
    if(!line.empty()) lines.push_back(std::move(line));
    return lines;
}

// Raw lines get blank space of the same width so tables line up with log messages
class PrefixFlag final : public spdlog::custom_flag_formatter {
  public:
    explicit PrefixFlag(bool timestamp) : m_timestamp{timestamp} {}

    void format(const spdlog::details::log_msg &msg, const std::tm &time,
                spdlog::memory_buf_t &dest) override {
        if(msg.logger_name == raw_logger_name) {
            Append(std::string(m_timestamp ? time_width + level_width : level_width, ' '), dest);
            return;
        }
        if(m_timestamp)
            fmt::format_to(std::back_inserter(dest), "{:02}:{:02}:{:02} ", time.tm_hour,
                           time.tm_min, time.tm_sec);
        msg.color_range_start = dest.size();
        Append(LevelName(msg.level), dest);
        msg.color_range_end = dest.size();
        dest.push_back(' ');
    }

    std::unique_ptr<custom_flag_formatter> clone() const override {
        return spdlog::details::make_unique<PrefixFlag>(m_timestamp);
    }

  private:
    bool m_timestamp;

    static void Append(std::string_view text, spdlog::memory_buf_t &dest) {
        dest.append(text.data(), text.data() + text.size());
    }

    static std::string_view LevelName(spdlog::level::level_enum level) {
        switch(level) {
        case spdlog::level::trace:
            return "TRACE";
        case spdlog::level::debug:
            return "DEBUG";
        case spdlog::level::info:
            return "INFO ";
        case spdlog::level::warn:
            return "WARN ";
        case spdlog::level::err:
            return "ERROR";
        case spdlog::level::critical:
            return "CRIT ";
        default:
            return "     ";
        }
    }
};

class ForwardingSink : public spdlog::sinks::sink {
  public:
    explicit ForwardingSink(std::shared_ptr<spdlog::sinks::sink> inner)
        : m_inner{std::move(inner)} {}
    void flush() override { m_inner->flush(); }
    void set_pattern(const std::string &pattern) override { m_inner->set_pattern(pattern); }
    void set_formatter(std::unique_ptr<spdlog::formatter> formatter) override {
        m_inner->set_formatter(std::move(formatter));
    }

  protected:
    std::shared_ptr<spdlog::sinks::sink> m_inner;
};

class LiveConsoleSink final : public ForwardingSink {
  public:
    using ForwardingSink::ForwardingSink;
    void log(const spdlog::details::log_msg &msg) override {
        auto &live = LiveLine::Instance();
        live.Suspend();
        m_inner->log(msg);
        live.Resume();
    }
};

// Keeps the log file plain text
class PlainTextSink final : public ForwardingSink {
  public:
    using ForwardingSink::ForwardingSink;
    void log(const spdlog::details::log_msg &msg) override {
        std::string_view payload(msg.payload.data(), msg.payload.size());
        if(payload.find('\033') == std::string_view::npos) {
            m_inner->log(msg);
            return;
        }
        const std::string plain = StripAnsi(payload);
        auto copy = msg;
        copy.payload = spdlog::string_view_t(plain.data(), plain.size());
        m_inner->log(copy);
    }
};

} // namespace

Style &achilles::console::GetStyle() {
    static Style style = DetectStyle();
    return style;
}

void achilles::console::Configure(bool batch_mode, bool force_ascii) {
    auto &style = GetStyle();
    if(batch_mode) style.interactive = false;
    if(force_ascii) style.unicode = false;
}

size_t achilles::console::DisplayWidth(std::string_view text) {
    const std::string plain = StripAnsi(text);
    return static_cast<size_t>(
        std::count_if(plain.begin(), plain.end(), [](char c) { return !IsContinuationByte(c); }));
}

size_t achilles::console::TerminalWidth() {
    winsize size{};
    if(ioctl(fileno(stdout), TIOCGWINSZ, &size) == 0 && size.ws_col > 0) return size.ws_col;
    if(const char *columns = std::getenv("COLUMNS")) {
        const long value = std::strtol(columns, nullptr, 10);
        if(value > 0) return static_cast<size_t>(value);
    }
    return 100;
}

std::string achilles::console::Pad(std::string_view text, size_t width, Align align) {
    const size_t length = DisplayWidth(text);
    if(length >= width) return std::string(text);
    const size_t fill = width - length;
    switch(align) {
    case Align::right:
        return std::string(fill, ' ') + std::string(text);
    case Align::center:
        return std::string(fill / 2, ' ') + std::string(text) + std::string(fill - fill / 2, ' ');
    case Align::left:
    default:
        return std::string(text) + std::string(fill, ' ');
    }
}

std::string achilles::console::Truncate(std::string_view text, size_t width) {
    if(DisplayWidth(text) <= width) return std::string(text);
    const std::string ellipsis = GetStyle().unicode ? "…" : "...";
    const size_t ellipsis_width = DisplayWidth(ellipsis);
    if(width <= ellipsis_width) return ellipsis.substr(0, width);

    // Escape sequences are dropped so a truncated color code cannot leak
    const std::string plain = StripAnsi(text);
    std::string result;
    size_t used = 0;
    for(size_t i = 0; i < plain.size(); ++i) {
        if(!IsContinuationByte(plain[i])) {
            if(used == width - ellipsis_width) break;
            ++used;
        }
        result += plain[i];
    }
    return result + ellipsis;
}

std::string achilles::console::FormatDuration(double seconds) {
    if(!std::isfinite(seconds) || seconds < 0) return "--";
    if(seconds < 1) return fmt::format("{:.0f}ms", seconds * 1000);
    if(seconds < 10) return fmt::format("{:.1f}s", seconds);
    auto total = static_cast<long long>(std::round(seconds));
    if(total < 60) return fmt::format("{}s", total);
    if(total < 3600) return fmt::format("{}m {:02}s", total / 60, total % 60);
    if(total < 86400) return fmt::format("{}h {:02}m", total / 3600, (total % 3600) / 60);
    return fmt::format("{}d {:02}h", total / 86400, (total % 86400) / 3600);
}

std::string achilles::console::FormatCount(double count) {
    const double magnitude = std::abs(count);
    if(magnitude < 1e3) return fmt::format("{:.0f}", count);
    if(magnitude < 1e6) return fmt::format("{:.3g}k", count / 1e3);
    if(magnitude < 1e9) return fmt::format("{:.3g}M", count / 1e6);
    return fmt::format("{:.3g}G", count / 1e9);
}

std::string achilles::console::FormatBytes(double bytes) {
    if(bytes < 1e3) return fmt::format("{:.0f} B", bytes);
    if(bytes < 1e6) return fmt::format("{:.3g} kB", bytes / 1e3);
    if(bytes < 1e9) return fmt::format("{:.3g} MB", bytes / 1e6);
    return fmt::format("{:.3g} GB", bytes / 1e9);
}

std::string achilles::console::FormatResult(double mean, double error) {
    return fmt::format("{:.5e} {} {:.2e}", mean, GetStyle().unicode ? "±" : "+/-", error);
}

std::vector<std::string> achilles::console::Wrap(const std::vector<std::string> &items,
                                                 size_t width, std::string_view separator) {
    std::vector<std::string> lines;
    std::string line;
    std::string_view line_end = separator;
    while(!line_end.empty() && line_end.back() == ' ') line_end.remove_suffix(1);

    for(size_t i = 0; i < items.size(); ++i) {
        const bool last = i + 1 == items.size();
        const std::string piece = items[i] + std::string(last ? "" : line_end);
        if(!line.empty() && DisplayWidth(line) + DisplayWidth(piece) > width) {
            while(!line.empty() && line.back() == ' ') line.pop_back();
            lines.push_back(line);
            line.clear();
        }
        line += items[i];
        if(!last) line += separator;
    }
    if(!line.empty()) lines.push_back(line);
    return lines;
}

std::string achilles::console::Colorize(std::string_view text, Color color) {
    if(!GetStyle().color || color == Color::none) return std::string(text);
    static constexpr int codes[] = {0, 1, 2, 31, 32, 33, 34, 35, 36};
    const int code = codes[static_cast<size_t>(color)];
    return fmt::format("\033[{}m{}\033[0m", code, text);
}

Table::Table(std::vector<std::string> headers, std::vector<Align> aligns)
    : m_headers{std::move(headers)}, m_aligns{std::move(aligns)} {
    m_aligns.resize(m_headers.size(), Align::left);
    m_min_widths.resize(m_headers.size(), 0);
}

void Table::SetMinWidths(std::vector<size_t> widths) {
    widths.resize(m_headers.size(), 0);
    m_min_widths = std::move(widths);
}

void Table::AddRow(std::vector<std::string> row) {
    row.resize(m_headers.size());
    m_rows.push_back(std::move(row));
}

void Table::AddSeparator() {
    m_rows.emplace_back();
}

size_t Table::MaxWidth() const {
    return m_max_width != 0 ? m_max_width : RawWidth();
}

size_t Table::InnerWidth(const std::vector<size_t> &widths) const {
    const size_t ncolumns = widths.size();
    // "│ a │ b │" for a grid, "│ a  b │" for a frame
    size_t inner = m_style == Style::grid ? 3 * ncolumns - 1 : 2 * ncolumns;
    for(const auto &width : widths) inner += width;
    return inner;
}

std::vector<size_t> Table::Widths() const {
    std::vector<size_t> widths(m_headers.size());
    for(size_t i = 0; i < m_headers.size(); ++i)
        widths[i] = std::max(m_min_widths[i], DisplayWidth(m_headers[i]));
    for(const auto &row : m_rows) {
        for(size_t i = 0; i < row.size(); ++i)
            widths[i] = std::max(widths[i], DisplayWidth(row[i]));
    }
    if(widths.empty()) return widths;

    if(!m_title.empty()) {
        const size_t needed = DisplayWidth(m_title) + 2;
        if(needed > InnerWidth(widths)) widths.back() += needed - InnerWidth(widths);
    }

    constexpr size_t min_width = 6;
    if(const size_t limit = MaxWidth(); limit != 0) {
        while(InnerWidth(widths) + 2 > limit) {
            auto widest = std::max_element(widths.begin(), widths.end());
            if(*widest <= min_width) break;
            --*widest;
        }
    }
    return widths;
}

std::string Table::RuleLine(const char *left, const char *mid, const char *right,
                            const std::vector<size_t> &widths) const {
    const auto &glyphs = TableGlyphs();
    std::string line = left;
    if(m_style == Style::frame) return line + Repeat(glyphs.horizontal, InnerWidth(widths)) + right;
    for(size_t i = 0; i < widths.size(); ++i) {
        line += Repeat(glyphs.horizontal, widths[i] + 2);
        line += i + 1 == widths.size() ? right : mid;
    }
    return line;
}

std::string Table::MiddleRule(const std::vector<size_t> &widths) const {
    const auto &glyphs = TableGlyphs();
    return RuleLine(glyphs.mid_left, glyphs.mid_mid, glyphs.mid_right, widths);
}

std::vector<std::string> Table::FormatRow(const std::vector<std::string> &row,
                                          const std::vector<size_t> &widths, bool header) const {
    if(!header && row.empty()) return {MiddleRule(widths)};
    const auto &glyphs = TableGlyphs();
    std::vector<std::vector<std::string>> cells;
    size_t nlines = 1;
    for(size_t i = 0; i < widths.size(); ++i) {
        cells.push_back(WrapText(i < row.size() ? row[i] : "", widths[i]));
        nlines = std::max(nlines, cells.back().size());
    }

    const std::string between =
        m_style == Style::grid ? std::string(" ") + glyphs.vertical + " " : std::string("  ");
    std::vector<std::string> lines;
    for(size_t line_index = 0; line_index < nlines; ++line_index) {
        std::string line = std::string(glyphs.vertical) + " ";
        for(size_t i = 0; i < widths.size(); ++i) {
            const std::string cell =
                line_index < cells[i].size() ? cells[i][line_index] : std::string();
            const Align align = header ? Align::left : m_aligns[i];
            std::string padded = Pad(cell, widths[i], align);
            if(header) padded = Colorize(padded, Color::bold);
            line += (i == 0 ? "" : between) + padded;
        }
        lines.push_back(line + " " + glyphs.vertical);
    }
    return lines;
}

std::vector<std::string> Table::HeaderLines() const {
    const auto widths = Widths();
    const auto &glyphs = TableGlyphs();
    std::vector<std::string> lines;
    if(m_title.empty()) {
        lines.push_back(RuleLine(glyphs.top_left, glyphs.top_mid, glyphs.top_right, widths));
    } else {
        const size_t inner = InnerWidth(widths);
        lines.push_back(RuleLine(glyphs.top_left, glyphs.horizontal, glyphs.top_right, widths));
        lines.push_back(std::string(glyphs.vertical) + " " +
                        Colorize(Pad(Truncate(m_title, inner - 2), inner - 2), Color::bold) + " " +
                        glyphs.vertical);
        lines.push_back(RuleLine(glyphs.mid_left, glyphs.top_mid, glyphs.mid_right, widths));
    }
    for(auto &line : FormatRow(m_headers, widths, true)) lines.push_back(std::move(line));
    lines.push_back(MiddleRule(widths));
    return lines;
}

std::string Table::RowLine(const std::vector<std::string> &row) const {
    std::string text;
    for(const auto &line : FormatRow(row, Widths(), false))
        text += (text.empty() ? "" : "\n") + line;
    return text;
}

std::string Table::FooterLine() const {
    const auto &glyphs = TableGlyphs();
    return RuleLine(glyphs.bottom_left, glyphs.bottom_mid, glyphs.bottom_right, Widths());
}

std::vector<std::string> Table::Render() const {
    auto lines = HeaderLines();
    const auto widths = Widths();
    for(const auto &row : m_rows) {
        for(auto &line : FormatRow(row, widths, false)) lines.push_back(std::move(line));
    }
    lines.push_back(FooterLine());
    return lines;
}

void Table::Print(Verbosity verbosity) const {
    achilles::console::Print(Render(), verbosity);
}

std::vector<std::string>
achilles::console::Box(const std::vector<std::vector<std::string>> &sections, size_t inner_width) {
    const auto &glyphs = BoxGlyphs();
    std::vector<std::string> lines;
    lines.push_back(std::string(glyphs.top_left) + Repeat(glyphs.horizontal, inner_width) +
                    glyphs.top_right);
    for(size_t i = 0; i < sections.size(); ++i) {
        if(i != 0)
            lines.push_back(std::string(glyphs.mid_left) + Repeat(glyphs.mid_mid, inner_width) +
                            glyphs.mid_right);
        for(const auto &line : sections[i])
            lines.push_back(std::string(glyphs.vertical) + Pad(line, inner_width) +
                            glyphs.vertical);
    }
    lines.push_back(std::string(glyphs.bottom_left) + Repeat(glyphs.horizontal, inner_width) +
                    glyphs.bottom_right);
    return lines;
}

std::vector<std::string>
achilles::console::Fields(const std::string &title,
                          const std::vector<std::pair<std::string, std::string>> &rows) {
    size_t width = 0;
    for(const auto &row : rows) width = std::max(width, DisplayWidth(row.first));
    std::vector<std::string> lines;
    if(!title.empty()) lines.push_back(Colorize(title, Color::bold));
    const std::string indent = title.empty() ? "" : "  ";
    // Long values (paths) wrap under themselves instead of back to the start of the line
    const size_t value_column = indent.size() + width + 2;
    const size_t limit = RawWidth();
    const size_t value_width = limit > value_column + 10 ? limit - value_column : 0;
    for(const auto &[name, value] : rows) {
        const auto parts = WrapText(value, value_width);
        lines.push_back(indent + Colorize(Pad(name, width), Color::dim) + "  " + parts.front());
        for(size_t i = 1; i < parts.size(); ++i)
            lines.push_back(std::string(value_column, ' ') + parts[i]);
    }
    return lines;
}

void achilles::console::Print(const std::vector<std::string> &lines, Verbosity verbosity) {
    // Raw lines have no level prefix, so "always" can use the warning level to get through
    // quiet mode without looking like a warning
    const auto level = verbosity == Verbosity::debug    ? spdlog::level::debug
                       : verbosity == Verbosity::always ? spdlog::level::warn
                                                        : spdlog::level::info;
    auto logger = spdlog::get(raw_logger_name);
    if(logger) {
        for(const auto &line : lines) logger->log(level, line);
        return;
    }
    // Without a logger there is no way to ask for debug output
    if(verbosity == Verbosity::debug) return;

    auto &live = LiveLine::Instance();
    live.Suspend();
    for(const auto &line : lines) fmt::print("{}\n", line);
    std::fflush(stdout);
    live.Resume();
}

void achilles::console::Print(const std::string &line, Verbosity verbosity) {
    Print(std::vector<std::string>{line}, verbosity);
}

LiveLine &LiveLine::Instance() {
    static LiveLine live;
    return live;
}

void LiveLine::Write(const std::string &text) {
    std::fputs(text.c_str(), stdout);
    std::fflush(stdout);
}

void LiveLine::Draw(const std::string &text) {
    std::lock_guard<std::recursive_mutex> lock(m_mutex);
    if(!GetStyle().interactive) return;
    m_text = text;
    Write("\r\033[2K" + m_text);
    m_drawn = true;
}

void LiveLine::Clear() {
    std::lock_guard<std::recursive_mutex> lock(m_mutex);
    Suspend();
    m_text.clear();
}

void LiveLine::Suspend() {
    std::lock_guard<std::recursive_mutex> lock(m_mutex);
    if(m_drawn) Write("\r\033[2K");
    m_drawn = false;
}

void LiveLine::Resume() {
    std::lock_guard<std::recursive_mutex> lock(m_mutex);
    if(m_text.empty() || !GetStyle().interactive) return;
    Write("\r\033[2K" + m_text);
    m_drawn = true;
}

namespace {
int terminal_level = spdlog::level::info;
} // namespace

int achilles::console::ConsoleLevel() {
    return terminal_level;
}

void achilles::console::SetupLoggers(int level, int log_level, int flush_time,
                                     const std::string &destination) {
    terminal_level = level;
    auto console_level = static_cast<spdlog::level::level_enum>(level);
    auto console_sink =
        std::make_shared<LiveConsoleSink>(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
    console_sink->set_level(console_level);

    auto file_level = static_cast<spdlog::level::level_enum>(log_level);
    auto file_sink = std::make_shared<PlainTextSink>(
        std::make_shared<spdlog::sinks::basic_file_sink_mt>(destination, true));
    file_sink->set_level(file_level);

    spdlog::sinks_init_list sink_list = {file_sink, console_sink};
    auto logger = std::make_shared<spdlog::logger>("achilles", sink_list);
    // The sinks filter for the terminal and the file, so the logger passes what either wants
    logger->set_level(std::min(console_level, file_level));
    logger->flush_on(spdlog::level::warn);

    // Shares the sinks, so tables end up in both the terminal and the log file
    auto raw_logger = std::make_shared<spdlog::logger>(raw_logger_name, sink_list);
    // The sinks decide what is shown, so debug tables can still reach a debug log file
    raw_logger->set_level(spdlog::level::trace);

    spdlog::drop("achilles");
    spdlog::drop(raw_logger_name);
    spdlog::register_logger(logger);
    spdlog::register_logger(raw_logger);

    // Both loggers share the sinks, so the formats are set on the sinks themselves
    auto make_formatter = [](bool timestamp) {
        auto formatter = spdlog::details::make_unique<spdlog::pattern_formatter>();
        formatter->add_flag<PrefixFlag>('*', timestamp).set_pattern("%*%v");
        // The time comes from a custom flag, which spdlog does not detect on its own
        formatter->need_localtime();
        return formatter;
    };
    file_sink->set_formatter(make_formatter(true));
    console_sink->set_formatter(make_formatter(false));
    spdlog::set_default_logger(logger);
    spdlog::flush_every(std::chrono::seconds(flush_time));
}
