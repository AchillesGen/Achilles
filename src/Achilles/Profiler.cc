#include "Achilles/Profiler.hh"
#include "Achilles/Logging.hh"
#include "fmt/chrono.h"
#include "fmt/format.h"

#include <algorithm>
#include <vector>

using achilles::Profiler;

Profiler &Profiler::Instance() {
    static Profiler profiler;
    return profiler;
}

void Profiler::Start() {
    m_start = Clock::now();
    spdlog::info("Start Time: {}", FormatTimestamp(std::chrono::system_clock::now()));
}

void Profiler::Record(const std::string &name, double seconds) {
    auto &section = m_sections[name];
    if(section.calls == 0) section.order = m_sections.size() - 1;
    section.time += seconds;
    section.calls++;
}

void Profiler::Reset() {
    m_sections.clear();
    m_start = Clock::now();
}

double Profiler::Time(const std::string &name) const {
    const auto it = m_sections.find(name);
    return it == m_sections.end() ? 0 : it->second.time;
}

size_t Profiler::Calls(const std::string &name) const {
    const auto it = m_sections.find(name);
    return it == m_sections.end() ? 0 : it->second.calls;
}

std::string Profiler::FormatTimestamp(std::chrono::system_clock::time_point time) {
    return fmt::format("{:%a %b %d %H:%M:%S %Y}",
                       fmt::localtime(std::chrono::system_clock::to_time_t(time)));
}

std::string Profiler::FormatTime(double seconds) {
    if(seconds < 1e-3) return fmt::format("{:.1f} us", seconds * 1e6);
    if(seconds < 1) return fmt::format("{:.1f} ms", seconds * 1e3);
    if(seconds < 60) return fmt::format("{:.3f} s", seconds);

    auto total = static_cast<int>(seconds);
    const int secs = total % 60;
    total /= 60;
    const int mins = total % 60;
    const int hours = total / 60;
    if(hours == 0) return fmt::format("{}m {:02d}s", mins, secs);
    return fmt::format("{}h {:02d}m {:02d}s", hours, mins, secs);
}

void Profiler::Summary() const {
    if(m_sections.empty()) return;

    // The reference point is the start of the run, so sections that are not timed
    // (setup, cascade, output, ...) show up in the unprofiled remainder
    const double total = std::chrono::duration<double>(Clock::now() - m_start).count();

    std::vector<const std::pair<const std::string, Section> *> ordered;
    ordered.reserve(m_sections.size());
    for(const auto &section : m_sections) ordered.push_back(&section);
    std::sort(ordered.begin(), ordered.end(), [](const auto *lhs, const auto *rhs) {
        return lhs->second.order < rhs->second.order;
    });

    const auto fraction = [&](double time) { return total > 0 ? 100 * time / total : 0.0; };

    spdlog::info("Timing summary:");
    spdlog::info("  {:<28} {:>6} {:>14} {:>8}", "Section", "Calls", "Time", "Fraction");
    double profiled = 0;
    for(const auto *section : ordered) {
        profiled += section->second.time;
        spdlog::info("  {:<28} {:>6} {:>14} {:>7.2f}%", section->first, section->second.calls,
                     FormatTime(section->second.time), fraction(section->second.time));
    }
    spdlog::info("  {:<28} {:>6} {:>14} {:>7.2f}%", "Other", "", FormatTime(total - profiled),
                 fraction(total - profiled));
    spdlog::info("  {:<28} {:>6} {:>14} {:>7.2f}%", "Total", "", FormatTime(total),
                 fraction(total));
    spdlog::info("End Time: {}", FormatTimestamp(std::chrono::system_clock::now()));
}
