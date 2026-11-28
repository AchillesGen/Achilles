// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/EventGen.hh"
#include "Achilles/Logging.hh"
#include "Achilles/Logo.hh"
#include "Achilles/Progress.hh"
#include "Achilles/ReferenceHandler.hh"
#include "Achilles/System.hh"
#include "Achilles/Version.hh"

#include "Plugins/Manager/PluginManager.hh"

#include <chrono>
#include <ctime>
#include <dlfcn.h>
#include <filesystem>

namespace fs = std::filesystem;
using namespace achilles::SystemVariables;

namespace achilles {

/** Gets the current time, logs it, and returns it as a number of seconds since epoch. */
static time_t logTime(std::string message) {
    time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    message += ctime(&time);
    // ctime ends with a newline, the logger adds its own
    if(!message.empty() && message.back() == '\n') message.pop_back();
    spdlog::info(message);
    return time;
}

void InitializePaths(const fs::path &root, const fs::path &libs, const fs::path &data,
                     const fs::path &flux) {
    achilles::PathVariables::Instance().InitializePython(root, libs, data, flux);
}

static bool initialized = false;
static std::string logFilePath;

void InitializeLogging(int verbosity, int log_verbosity, std::string logfile) {
    // If a python user calls initialize_logging() more than once, they should be warned
    if(initialized) {
        spdlog::warn("Logger has already been initialized, ignoring this call");
        return;
    }

    logFilePath = logfile;
    CreateLogger(verbosity, log_verbosity, 1, logfile);
    GitInformation();

    // Install signal handlers
    std::signal(SIGTERM, SignalHandler);
    std::signal(SIGSEGV, SignalHandler);
    std::signal(SIGINT, SignalHandler);
    std::signal(SIGABRT, SignalHandler);

    // Ensure reference handle is initialized
    achilles::Reference main_ref{achilles::ReferenceType::inspire, "Isaacson:2022cwh",
                                 "Main Achilles reference"};
    achilles::ReferenceHandler &ref_handler = achilles::ReferenceHandler::Handle();
    ref_handler.AddReference(main_ref);

    achilles::Plugin::Manager plugin_manager;
    initialized = true;
}

void GenerateEvents(const std::string &runcard, std::vector<std::string> &shargs, bool batchMode) {
    // Python users should be able to call "generate_events()" without needing to initialize_logging
    // first.
    if(!initialized) InitializeLogging(2, 2, "achilles.log");

    if(!fs::exists(runcard)) {
        if(runcard == "run.yml") {
            spdlog::error("Achilles: Could not find \"run.yml\". Copying over default run card to "
                          "this location");
            // Resolve the default run card through the standard search path
            try {
                auto default_runcard =
                    achilles::Filesystem::FindFile("data/default/run.yml", "Achilles");
                fs::copy(default_runcard, fs::current_path() / "run.yml");
            } catch(const std::exception &error) {
                spdlog::error("Achilles: Could not copy the default run card: {}", error.what());
            }
        } else
            spdlog::error("Achilles: Could not find \"" + runcard + "\".");
        return;
    }

    logTime("Start Time: ");
    Stopwatch run_clock;
    console::Configure(batchMode);
    auto &tracker = RunTracker::Instance();
    tracker.Reset();

    std::string success = "Failed";
    try {
        tracker.AddStep("setup", "Setup");
        tracker.SetStepState("setup", StepState::running);
        tracker.BeginStage("Setup");
        achilles::EventGen generator(runcard, shargs);
        tracker.EndStage("Read " + runcard);
        tracker.SetStepState("setup", StepState::done);
        generator.Initialize();
        generator.GenerateEvents(batchMode);
        success = "Success";
    } catch(const std::runtime_error &error) {
        spdlog::error(error.what());
        tracker.EndStage("Stopped by an error");
    }

    spdlog::debug("Event Run Concluded - " + success);

    // The event file is only complete once the generator is destroyed
    if(!tracker.OutputFile().empty()) {
        std::string output = tracker.OutputFile();
        if(!fs::exists(output) && fs::exists(output + ".gz")) output += ".gz";
        tracker.SetResult(
            "Output",
            fs::exists(output)
                ? fmt::format("{} ({})", output,
                              console::FormatBytes(static_cast<double>(fs::file_size(output))))
                : output + " (not written)");
    }
    tracker.SetResult("Log file", logFilePath);
    tracker.SetResult("References", "References.txt (papers to cite for this run)");
    tracker.PrintSummary(run_clock.Elapsed(), success);
}

} // namespace achilles
