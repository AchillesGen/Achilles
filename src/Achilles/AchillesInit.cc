#include "Achilles/Achilles.hh"
#include "Achilles/EventGen.hh"
#include "Achilles/Logging.hh"
#include "Achilles/Logo.hh"
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
    spdlog::info(message);
    return time;
}

/** Puts a potentially-large number of seconds into a more human-readable form
 *  There might've been a library for this but I coded it myself anyway -Hayden */
static std::string formatTime(time_t seconds) {
    std::string output = std::to_string(seconds % 60) + "s";
    seconds /= 60;
    if(seconds == 0) return output;
    output = std::to_string(seconds % 60) + "m " + output;
    seconds /= 60;
    if(seconds == 0) return output;
    output = std::to_string(seconds % 24) + "h " + output;
    seconds /= 24;
    if(seconds == 0) return output;
    return std::to_string(seconds) + "d " + output;
}

void InitializePaths(const fs::path &root, const fs::path &libs, const fs::path &data,
                     const fs::path &flux) {
    achilles::PathVariables::Instance().InitializePython(root, libs, data, flux);
}

void GenerateEvents(const std::string &runcard, std::vector<std::string> &shargs, int verbosity,
                    int log_verbosity, const std::string &logFilePath, bool batchMode) {
    CreateLogger(verbosity, log_verbosity, 1, logFilePath);
    GitInformation();

    // Install signal handlers
    std::signal(SIGTERM, SignalHandler);
    std::signal(SIGSEGV, SignalHandler);
    std::signal(SIGINT, SignalHandler);
    std::signal(SIGABRT, SignalHandler);

    // Ensure reference handle is initialized
    achilles::Reference main_ref{achilles::ReferenceType::inspire, "Isaacson:2022cwh",
                                 "Main Achilles reference"};
    auto &ref_handler = achilles::ReferenceHandler::Handle();
    ref_handler.AddReference(main_ref);

    if(!fs::exists(runcard)) {
        if(runcard == "run.yml") {
            spdlog::error("Achilles: Could not find \"run.yml\". Copying over default run card to "
                          "this location");
            fs::copy(achilles::PathVariables::Instance().DataPath() / "default/run.yml",
                     fs::current_path());
        } else
            spdlog::error("Achilles: Could not find \"" + runcard + "\".");
        return;
    }

    achilles::Plugin::Manager plugin_manager;

    time_t startTime = logTime("Start Time: ");

    std::string success = "Failed.";
    try {
        achilles::EventGen generator(runcard, shargs);
        generator.Initialize();
        generator.GenerateEvents(batchMode);
        success = "Success!";
    } catch(const std::runtime_error &error) { spdlog::error(error.what()); }

    spdlog::info("Event Run Concluded - " + success);
    spdlog::info("Records of this run can be found in \"" + logFilePath + "\"");

    time_t endTime = logTime("End Time: ");
    spdlog::info("Run Duration: " + formatTime(endTime - startTime));
}

} // namespace achilles
