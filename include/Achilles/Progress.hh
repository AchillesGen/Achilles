// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef PROGRESS_HH
#define PROGRESS_HH

#include <chrono>
#include <map>
#include <string>
#include <vector>

namespace achilles {

class Stopwatch {
  public:
    Stopwatch() : m_start{std::chrono::steady_clock::now()} {}
    void Reset() { m_start = std::chrono::steady_clock::now(); }
    double Elapsed() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - m_start).count();
    }

  private:
    std::chrono::steady_clock::time_point m_start;
};

// Redrawn in place on a terminal; otherwise logs every 10% or every minute
//   ● Setup  ● G1  G2: ◉ Optimize ━━━━━━╸───── 52%  ○ Max weight  ○ Events  ETA ~3m 41s
class ProgressBar {
  public:
    ProgressBar(std::string label, double total);
    ProgressBar(const ProgressBar &) = delete;
    ProgressBar &operator=(const ProgressBar &) = delete;
    ~ProgressBar();

    void SetTotal(double total, bool estimated = false);
    // Off when the stage reports progress another way when not drawn in place
    void SetLogLines(bool log) { m_log_lines = log; }
    void Update(double done, std::string detail = "");
    void Finish();

    double Done() const { return m_done; }
    double Total() const { return m_total; }
    double Elapsed() const { return m_clock.Elapsed(); }
    // Seconds left in this stage, negative when it cannot be estimated yet
    double Remaining() const;

    std::string Status() const;

  private:
    void Draw();
    bool DrawSteps(size_t width);
    void Log();
    std::string Percent() const;
    std::string Timing() const;
    std::string EtaText() const;

    std::string m_label, m_detail;
    double m_total{}, m_done{};
    bool m_estimated{false}, m_finished{false}, m_log_lines{true};
    Stopwatch m_clock;
    double m_last_draw{-1}, m_last_log{0}, m_last_log_fraction{0};
    // Layout only gets more compact at a fixed width, so the bar does not jump
    size_t m_width{0}, m_layout{0}, m_bar_width{0};
};

struct StageRecord {
    std::string name;
    double seconds;
    std::string detail;
    bool cached{false};
};

enum class StepState { pending, running, done, cached };

struct RunStep {
    std::string key, label, group;
    StepState state{StepState::pending};
};

struct CrossSectionRow {
    std::string name;
    double xsec, error;
    double share, efficiency;
};

struct TimeEstimate {
    double seconds{0};
    // False when seconds is only a lower bound
    bool complete{true};
};

// Pending work is counted in integrand calls and converted with the last measured cost per call
class RunTracker {
  public:
    static RunTracker &Instance();

    void Reset();
    void BeginStage(std::string name);
    void EndStage(std::string detail = "");
    void AddCachedStage(std::string name, std::string detail);

    // A negative number of calls marks work whose size is not known yet
    void SetPending(const std::string &key, double calls);
    void RemovePending(const std::string &key);
    void SetSecondsPerCall(double seconds) { m_seconds_per_call = seconds; }
    double SecondsPerCall() const { return m_seconds_per_call; }
    TimeEstimate Pending() const;

    // Adding a key twice has no effect
    void AddStep(const std::string &key, std::string label, std::string group = "");
    void SetStepState(const std::string &key, StepState state);
    const std::vector<RunStep> &Steps() const { return m_steps; }

    const std::vector<StageRecord> &Stages() const { return m_stages; }

    void AddCrossSection(CrossSectionRow row) { m_cross_sections.push_back(std::move(row)); }
    // Setting a name again replaces its value in place
    void SetResult(const std::string &name, std::string value);
    const std::vector<std::pair<std::string, std::string>> &Results() const { return m_results; }
    void SetOutputFile(std::string path) { m_output_file = std::move(path); }
    const std::string &OutputFile() const { return m_output_file; }

    void PrintSummary(double total_seconds, const std::string &status) const;

  private:
    RunTracker() = default;

    std::vector<StageRecord> m_stages;
    std::string m_current;
    Stopwatch m_stage_clock;
    std::map<std::string, double> m_pending;
    std::vector<RunStep> m_steps;
    std::vector<CrossSectionRow> m_cross_sections;
    std::vector<std::pair<std::string, std::string>> m_results;
    std::string m_output_file;
    double m_seconds_per_call{-1};
};

} // namespace achilles

#endif
