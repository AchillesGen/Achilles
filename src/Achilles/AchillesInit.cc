#include "Achilles/Achilles.hh"
#include "Achilles/EventGen.hh"
#include "Achilles/Logging.hh"
#include "Achilles/Logo.hh"
#include "Achilles/Profiler.hh"
#include "Achilles/ReferenceHandler.hh"
#include "Achilles/System.hh"
#include "Achilles/Version.hh"

#include "Plugins/Manager/PluginManager.hh"

#include <dlfcn.h>
#include <filesystem>

namespace fs = std::filesystem;
using namespace achilles::SystemVariables;

namespace achilles {

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
            fs::copy(achilles::PathVariables::Instance().DataPath() / "default/run.yml",
                     fs::current_path());
        } else
            spdlog::error("Achilles: Could not find \"" + runcard + "\".");
        return;
    }

    Profiler::Instance().Start();

    std::string success = "Failed.";
    try {
        achilles::EventGen generator(runcard, shargs);
        generator.Initialize();
        generator.GenerateEvents(batchMode);
        success = "Success!";
    } catch(const std::runtime_error &error) { spdlog::error(error.what()); }

    Profiler::Instance().Summary();

    spdlog::info("Event Run Concluded - " + success);
    spdlog::info("Records of this run can be found in \"" + logFilePath + "\"");
}

} // namespace achilles
