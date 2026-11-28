// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Progress.hh"
#include "Achilles/Console.hh"

#include "fmt/format.h"
#include "fmt/ranges.h"
#include "spdlog/spdlog.h"

#include <algorithm>
#include <cmath>
#include <iterator>

using namespace achilles;
using achilles::console::Color;
using achilles::console::Colorize;
using achilles::console::FormatDuration;

namespace {

constexpr size_t label_width = 12;
constexpr double redraw_interval = 0.1;
constexpr double log_interval = 60;
constexpr double log_fraction = 0.1;

bool DrawLive() {
    // Debug output scrolls too quickly for a live line to be useful
    return console::GetStyle().interactive && console::ConsoleLevel() >= spdlog::level::info;
}

std::string Separator() {
    return console::GetStyle().unicode ? " · " : " | ";
}

std::string BarGlyphs(size_t width, double fraction) {
    fraction = std::clamp(fraction, 0.0, 1.0);
    const auto filled = static_cast<size_t>(fraction * static_cast<double>(width));
    const bool head = filled < width && fraction > 0;
    const size_t empty = width - filled - (head ? 1 : 0);

    std::string done, rest;
    if(console::GetStyle().unicode) {
        for(size_t i = 0; i < filled; ++i) done += "━";
        if(head) done += "╸";
        for(size_t i = 0; i < empty; ++i) rest += "─";
        return Colorize(done, Color::magenta) + Colorize(rest, Color::dim);
    }
    done = std::string(filled, '=') + (head ? ">" : "");
    return "[" + done + std::string(empty, ' ') + "]";
}

} // namespace

ProgressBar::ProgressBar(std::string label, double total)
    : m_label{std::move(label)}, m_total{total} {}

ProgressBar::~ProgressBar() {
    Finish();
}

void ProgressBar::SetTotal(double total, bool estimated) {
    m_total = total;
    m_estimated = estimated;
}

void ProgressBar::Update(double done, std::string detail) {
    if(m_finished) return;
    m_done = done;
    m_detail = std::move(detail);

    const double now = Elapsed();
    if(DrawLive()) {
        if(now - m_last_draw < redraw_interval) return;
        m_last_draw = now;
        Draw();
        return;
    }

    if(!m_log_lines) return;
    const double fraction = m_total > 0 ? m_done / m_total : 0;
    if(m_done > 0 &&
       (fraction - m_last_log_fraction >= log_fraction || now - m_last_log >= log_interval)) {
        m_last_log = now;
        m_last_log_fraction = fraction;
        Log();
    }
}

void ProgressBar::Finish() {
    if(m_finished) return;
    m_finished = true;
    console::LiveLine::Instance().Clear();
}

double ProgressBar::Remaining() const {
    const double elapsed = Elapsed();
    if(m_done <= 0 || elapsed <= 0) return -1;
    const double rate = m_done / elapsed;
    return std::max(m_total - m_done, 0.0) / rate;
}

std::string ProgressBar::Percent() const {
    if(m_total <= 0) return "";
    return fmt::format("{:>3.0f}%", 100 * std::clamp(m_done / m_total, 0.0, 1.0));
}

std::string ProgressBar::Status() const {
    std::string status = Percent();
    if(!m_detail.empty()) status += "  " + m_detail;
    return status + Separator() + Timing();
}

std::string ProgressBar::Timing() const {
    return FormatDuration(Elapsed()) + Separator() + EtaText();
}

std::string ProgressBar::EtaText() const {
    const std::string sep = Separator();
    std::string status;

    const double remaining = Remaining();
    const std::string approx = m_estimated ? "~" : "";
    status += "ETA " + (remaining < 0 ? "--" : approx + FormatDuration(remaining));

    const auto pending = RunTracker::Instance().Pending();
    if(remaining >= 0 && (pending.seconds > 0 || !pending.complete)) {
        std::string bound = pending.complete ? "~" : (console::GetStyle().unicode ? "≥" : ">=");
        status += sep + "run " + bound + FormatDuration(remaining + pending.seconds);
    }
    return status;
}

void ProgressBar::Draw() {
    const size_t width = console::TerminalWidth();
    if(width != m_width) {
        m_width = width;
        m_layout = 0;
        m_bar_width = 0;
    }
    if(DrawSteps(width)) return;

    const std::string label = Colorize(console::Pad(m_label, label_width), Color::bold);
    const std::string status = Status();

    const size_t used = label_width + 2 + console::DisplayWidth(status) + 1;
    size_t bar_width = used + 10 < width ? std::min<size_t>(width - used, 30) : 10;
    if(m_bar_width != 0) bar_width = std::min(bar_width, m_bar_width);
    m_bar_width = bar_width;
    const double fraction = m_total > 0 ? m_done / m_total : 0;

    std::string line = label + " " + BarGlyphs(bar_width, fraction) + " " + status;
    if(console::DisplayWidth(line) >= width) line = console::Truncate(line, width - 1);
    console::LiveLine::Instance().Draw(line);
}

bool ProgressBar::DrawSteps(size_t width) {
    const auto &steps = RunTracker::Instance().Steps();
    const bool running = std::any_of(steps.begin(), steps.end(), [](const RunStep &step) {
        return step.state == StepState::running;
    });
    if(steps.size() < 2 || !running) return false;

    std::vector<std::vector<const RunStep *>> items;
    for(const auto &step : steps) {
        if(!items.empty() && !step.group.empty() && items.back().front()->group == step.group)
            items.back().push_back(&step);
        else
            items.push_back({&step});
    }

    const bool unicode = console::GetStyle().unicode;
    const double fraction = m_total > 0 ? m_done / m_total : 0;
    const std::string sep = Separator();
    const std::string eta = EtaText();
    const std::string elapsed = FormatDuration(Elapsed());
    auto dial = [&](StepState state) {
        switch(state) {
        case StepState::done:
            return Colorize(unicode ? "●" : "[x]", Color::green);
        case StepState::running:
            return Colorize(unicode ? "◉" : "[>]", Color::magenta);
        case StepState::cached:
            return Colorize(unicode ? "◌" : "[c]", Color::cyan);
        case StepState::pending:
        default:
            return Colorize(unicode ? "○" : "[ ]", Color::dim);
        }
    };

    // Most to least descriptive; pending labels are kept longest
    struct Layout {
        bool done_labels, pending_labels, elapsed, detail;
        size_t bar_width;
    };
    static constexpr Layout layouts[] = {
        {true, true, true, true, 20},    {false, true, true, true, 20},
        {false, true, false, true, 16},  {false, true, false, false, 12},
        {false, false, false, true, 12}, {false, false, false, false, 12},
    };

    std::string line;
    for(; m_layout < std::size(layouts); ++m_layout) {
        const Layout &layout = layouts[m_layout];
        // "bare": only the dial is shown
        auto step_text = [&](const RunStep &step, bool &bare) {
            bare = false;
            switch(step.state) {
            case StepState::running:
                return dial(step.state) + " " + Colorize(step.label, Color::bold) + " " +
                       BarGlyphs(layout.bar_width, fraction) + " " + Percent();
            case StepState::done:
                bare = !layout.done_labels;
                return dial(step.state) + (bare ? "" : " " + step.label);
            case StepState::cached:
                bare = !layout.done_labels;
                return dial(step.state) +
                       (bare ? "" : " " + step.label + Colorize(" (cached)", Color::dim));
            case StepState::pending:
            default:
                bare = !layout.pending_labels;
                return dial(step.state) + (bare ? "" : " " + Colorize(step.label, Color::dim));
            }
        };
        // Single space after a bare dial in a group, else only between consecutive bare dials
        auto join = [](const std::vector<std::pair<std::string, bool>> &parts, bool tight = false) {
            std::string result;
            for(size_t i = 0; i < parts.size(); ++i) {
                const bool single = parts[i - (i != 0)].second && (tight || parts[i].second);
                if(i != 0) result += single ? " " : "  ";
                result += parts[i].first;
            }
            return result;
        };

        std::vector<std::pair<std::string, bool>> parts;
        for(const auto &item : items) {
            const RunStep &first = *item.front();
            bool bare = false;
            if(first.group.empty()) {
                const std::string text = step_text(first, bare);
                parts.emplace_back(text, bare);
                continue;
            }

            auto all = [&](auto predicate) {
                return std::all_of(item.begin(), item.end(), predicate);
            };
            const bool group_running =
                !all([](const RunStep *step) { return step->state != StepState::running; });
            const bool all_cached =
                all([](const RunStep *step) { return step->state == StepState::cached; });
            const bool all_finished = all([](const RunStep *step) {
                return step->state == StepState::done || step->state == StepState::cached;
            });
            // Only the running group shows its steps
            if(!group_running && all_cached) {
                parts.emplace_back(dial(StepState::cached) + " " + first.group +
                                       Colorize(" cached", Color::dim),
                                   false);
            } else if(!group_running && all_finished) {
                parts.emplace_back(dial(StepState::done) + " " + first.group, false);
            } else if(!group_running) {
                parts.emplace_back(
                    dial(StepState::pending) + " " + Colorize(first.group, Color::dim), false);
            } else {
                std::vector<std::pair<std::string, bool>> sub;
                for(const auto *step : item) {
                    const std::string text = step_text(*step, bare);
                    sub.emplace_back(text, bare);
                }
                parts.emplace_back(Colorize(first.group + ":", Color::bold) + " " + join(sub, true),
                                   false);
            }
        }
        line = join(parts);

        std::vector<std::string> tail;
        if(layout.detail && !m_detail.empty()) tail.push_back(m_detail);
        if(layout.elapsed) tail.push_back(elapsed);
        tail.push_back(eta);
        line += fmt::format("   {}", fmt::join(tail, sep));
        if(console::DisplayWidth(line) < width) break;
    }
    m_layout = std::min(m_layout, std::size(layouts) - 1);

    if(console::DisplayWidth(line) >= width) line = console::Truncate(line, width - 1);
    console::LiveLine::Instance().Draw(line);
    return true;
}

void ProgressBar::Log() {
    spdlog::info("{}: {}", m_label, Status());
}

RunTracker &RunTracker::Instance() {
    static RunTracker tracker;
    return tracker;
}

void RunTracker::Reset() {
    m_stages.clear();
    m_steps.clear();
    m_cross_sections.clear();
    m_results.clear();
    m_output_file.clear();
    m_current.clear();
    m_pending.clear();
    m_seconds_per_call = -1;
}

void RunTracker::BeginStage(std::string name) {
    if(!m_current.empty()) EndStage();
    m_current = std::move(name);
    m_stage_clock.Reset();
}

void RunTracker::EndStage(std::string detail) {
    if(m_current.empty()) return;
    m_stages.push_back({m_current, m_stage_clock.Elapsed(), std::move(detail)});
    m_current.clear();
}

void RunTracker::AddStep(const std::string &key, std::string label, std::string group) {
    for(const auto &step : m_steps)
        if(step.key == key) return;
    m_steps.push_back({key, std::move(label), std::move(group), StepState::pending});
}

void RunTracker::SetStepState(const std::string &key, StepState state) {
    for(auto &step : m_steps)
        if(step.key == key) step.state = state;
}

void RunTracker::SetPending(const std::string &key, double calls) {
    m_pending[key] = calls;
}

void RunTracker::RemovePending(const std::string &key) {
    m_pending.erase(key);
}

TimeEstimate RunTracker::Pending() const {
    TimeEstimate estimate;
    for(const auto &[key, calls] : m_pending) {
        if(calls == 0) continue;
        if(calls < 0 || m_seconds_per_call < 0) {
            estimate.complete = false;
            continue;
        }
        estimate.seconds += calls * m_seconds_per_call;
    }
    return estimate;
}

void RunTracker::AddCachedStage(std::string name, std::string detail) {
    m_stages.push_back({std::move(name), 0, std::move(detail), true});
}

void RunTracker::SetResult(const std::string &name, std::string value) {
    for(auto &result : m_results) {
        if(result.first == name) {
            result.second = std::move(value);
            return;
        }
    }
    m_results.emplace_back(name, std::move(value));
}

void RunTracker::PrintSummary(double total_seconds, const std::string &status) const {
    // Shown even in quiet mode
    using console::Align;
    using console::Color;
    console::Table table({"Stage", "Time", "Details"}, {Align::left, Align::right, Align::left});
    table.SetStyle(console::Table::Style::frame);
    for(const auto &stage : m_stages)
        table.AddRow(
            {stage.name, stage.cached ? "cached" : FormatDuration(stage.seconds), stage.detail});
    table.AddSeparator();
    table.AddRow({"Total", FormatDuration(total_seconds), status});
    table.Print(console::Verbosity::always);

    if(!m_cross_sections.empty()) {
        console::Table xsecs(
            {"Process group", "XS [nb]", "AbsErr [nb]", "RelErr", "Events", "Unw. eff"},
            {Align::left, Align::right, Align::right, Align::right, Align::right, Align::right});
        xsecs.SetStyle(console::Table::Style::frame);
        for(const auto &row : m_cross_sections) {
            if(row.name == "Total") xsecs.AddSeparator();
            xsecs.AddRow({Colorize(row.name, Color::bold),
                          Colorize(fmt::format("{:.5e}", row.xsec), Color::blue),
                          Colorize(fmt::format("{:.2e}", row.error), Color::red),
                          Colorize(fmt::format("{:.2f} %", 100 * row.error / std::abs(row.xsec)),
                                   Color::red),
                          fmt::format("{:.1f} %", 100 * row.share),
                          fmt::format("{:.2f} %", 100 * row.efficiency)});
        }
        xsecs.Print(console::Verbosity::always);
    }

    if(!m_results.empty())
        console::Print(console::Fields("", m_results), console::Verbosity::always);
}
