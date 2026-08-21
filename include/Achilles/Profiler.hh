#ifndef ACHILLES_PROFILER_HH
#define ACHILLES_PROFILER_HH

#include <chrono>
#include <map>
#include <string>

namespace achilles {

/// Names of the sections timed during a standard event generation run
namespace ProfileSection {
constexpr auto integrator_optimization = "Integrator Optimization";
constexpr auto max_weight_scan = "Max Weight Scan";
constexpr auto event_generation = "Event Generation";
} // namespace ProfileSection

/// Accumulates the wall-clock time spent in named sections of a run
class Profiler {
  public:
    using Clock = std::chrono::steady_clock;

    struct Section {
        double time{};
        size_t calls{};
        size_t order{};
    };

    static Profiler &Instance();

    /// Mark the start of the run, both as a timestamp and as the reference point used
    /// for the fraction of the run spent in each section
    void Start();
    void Record(const std::string &, double);
    void Summary() const;
    /// Drop all recorded sections and restart the run clocks without logging
    void Reset();

    double Time(const std::string &) const;
    size_t Calls(const std::string &) const;
    const std::map<std::string, Section> &Sections() const { return m_sections; }

    static std::string FormatTime(double);
    static std::string FormatTimestamp(std::chrono::system_clock::time_point);

  private:
    Profiler() = default;

    Clock::time_point m_start{Clock::now()};
    std::map<std::string, Section> m_sections;
};

/// Times the enclosing scope and adds the result to the named Profiler section
class ScopedTimer {
  public:
    explicit ScopedTimer(std::string name) : m_name{std::move(name)} {}
    ~ScopedTimer() { Stop(); }
    ScopedTimer(const ScopedTimer &) = delete;
    ScopedTimer(ScopedTimer &&) = delete;
    ScopedTimer &operator=(const ScopedTimer &) = delete;
    ScopedTimer &operator=(ScopedTimer &&) = delete;

    double Elapsed() const {
        return std::chrono::duration<double>(Profiler::Clock::now() - m_start).count();
    }
    /// Record the time so far, and prevent the destructor from recording it again
    void Stop() {
        if(m_stopped) return;
        m_stopped = true;
        Profiler::Instance().Record(m_name, Elapsed());
    }

  private:
    std::string m_name;
    Profiler::Clock::time_point m_start{Profiler::Clock::now()};
    bool m_stopped{false};
};

} // namespace achilles

#endif // ACHILLES_PROFILER_HH
