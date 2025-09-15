#include "Achilles/EventGen.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Logging.hh"
#include "Achilles/Logo.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ReferenceHandler.hh"
#include "Achilles/System.hh"
#include "Achilles/Version.hh"
#include "Achilles/fortran/FNuclearModel.hh"
#include "Plugins/Manager/PluginManager.hh"
#include "git.h"
#ifdef ACHILLES_SHERPA_INTERFACE
#include "Plugins/Sherpa/Channels.hh"
#include "Plugins/Sherpa/SherpaInterface.hh"
#endif

#include "docopt.h"

#include <chrono>
#include <ctime>
#include <dlfcn.h>
#include <filesystem>

using namespace achilles::SystemVariables;
using namespace achilles::PathVariables;
namespace fs = std::filesystem;

void GitInformation() {
    std::string msg;
    if(git::IsPopulated()) {
        spdlog::info("Achilles git information");
        spdlog::info("    Commit: {}", git::CommitSHA1());
        spdlog::info("    Branch: {}", git::Branch());
        spdlog::info("    Local Changes: {}", git::AnyUncommittedChanges() ? "Yes" : "No");
    } else {
        spdlog::warn("This is not a git repository version of Achilles");
    }
}

static const std::string USAGE =
    R"(
    Usage:
      achilles [<input>] [-v | -vv] [-l | -ll] [-b] [--logfile=<logfile_name>] [-s | --sherpa=<sherpa>...]
      achilles --display-cuts
      achilles --display-ps
      achilles --display-ff
      achilles --display-int-models
      achilles --display-nuc-models
      achilles (-h | --help)
      achilles --version

    Options:
      -v[v]                                 Increase verbosity level.
      -l[l]                                 Increase log verbosity
                                              (Note: Log verbosity is never lower than total level)
      -b                                    Batch Mode (makes output more log-friendly)
      -h --help                             Show this screen.
      --version                             Show version.
      --logfile=<logfile_name>              Change the logging output destination
      -s <sherpa> --sherpa=<sherpa>         Define Sherpa option.
      --display-cuts                        Display the available cuts
      --display-ps                          Display the available phase spaces
      --display-ff                          Display the available form factors
      --display-int-models                  Display the available cascade interaction models
      --display-nuc-models                  Display the available nuclear interaction models
)";

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

void GenerateEvents(const std::string &runcard, const std::vector<std::string> &shargs,
                    bool batchMode) {
    achilles::EventGen generator(runcard, shargs);
    generator.Initialize();
    generator.GenerateEvents(batchMode);
}

int main(int argc, char *argv[]) {
    Splash();
    std::map<std::string, docopt::value> args =
        docopt::docopt(USAGE, {argv + 1, argv + argc},
                       true,                                          // show help if requested
                       fmt::format("achilles {}", ACHILLES_VERSION)); // version string

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

    auto verbosity = static_cast<int>(2 - args["-v"].asLong());
    auto log_verbosity = std::min(verbosity, static_cast<int>(2 - args["-l"].asLong()));

    bool batchMode = args["-b"].asBool();

    std::string logFilePath = "achilles.log";
    if(args["--logfile"].isString()) logFilePath = args["--logfile"].asString();

    CreateLogger(verbosity, log_verbosity, 1, logFilePath);
    GitInformation();
    achilles::Plugin::Manager plugin_manager;

    // Install signal handlers
    std::signal(SIGTERM, SignalHandler);
    std::signal(SIGSEGV, SignalHandler);
    std::signal(SIGINT, SignalHandler);
    std::signal(SIGABRT, SignalHandler);

    auto verbosity = static_cast<int>(2 - args["-v"].asLong());
    auto log_verbosity = std::min(verbosity, static_cast<int>(2 - args["-l"].asLong()));
    CreateLogger(verbosity, log_verbosity, 1);
    GitInformation();

    if(args["--display-cuts"].asBool()) {
        achilles::CutFactory<achilles::OneParticleCut>::Display();
        achilles::CutFactory<achilles::TwoParticleCut>::Display();
        return 0;
    }

    if(args["--display-ps"].asBool()) {
        achilles::Factory<achilles::HadronicBeamMapper, const achilles::ProcessInfo &,
                          size_t>::Display();
        achilles::Factory<achilles::FinalStateMapper, std::vector<double>>::Display();
#ifdef ACHILLES_SHERPA_INTERFACE
        achilles::Factory<PHASIC::Channels, std::vector<double>>::Display();
#endif // ACHILLES_SHERPA_INTERFACE
        return 0;
    }

    if(args["--display-ff"].asBool()) {
        achilles::FormFactorFactory::Display();
        return 0;
    }

    if(args["--display-int-models"].asBool()) {
        achilles::InteractionFactory::Display();
        return 0;
    }

    achilles::FortranModel::RegisterModels();
    if(args["--display-nuc-models"].asBool()) {
        achilles::NuclearModelFactory::Display();
        achilles::FortranModel::DisplayModels();
        return 0;
    }

    std::string runcard = "run.yml";
    if(args["<input>"].isString()) {
        runcard = args["<input>"].asString();
        if(!fs::exists(runcard)) {
            spdlog::error("Achilles: Could not find \"" + runcard + "\".");
            return 1;
        }
    } else {
        // Ensure file exists, otherwise copy template file to current location
        if(!fs::exists(runcard)) {
            spdlog::error("Achilles: Could not find \"run.yml\". Copying over default run card to "
                          "this location");
            if(!fs::exists(achilles::PathVariables::installData)) {
                fs::copy(achilles::PathVariables::buildData / fs::path("default/run.yml"),
                         fs::current_path());
            } else {
                fs::copy(achilles::PathVariables::installData / fs::path("default/run.yml"),
                         fs::current_path());
            }
            return 1;
        }
    }

    std::vector<std::string> shargs;
    if(args["--sherpa"].isStringList()) shargs = args["--sherpa"].asStringList();

    time_t startTime = logTime("Start Time: ");

    std::string success = "Failed.";
    try {
        GenerateEvents(runcard, shargs, batchMode);
        success = "Success!";
    } catch(const std::runtime_error &error) { spdlog::error(error.what()); }
    spdlog::info("Event Run Concluded - " + success);
    spdlog::info("Records of this run can be found in \"" + logFilePath + "\"");

    time_t endTime = logTime("End Time: ");
    spdlog::info("Run Duration: " + formatTime(endTime - startTime));

    return 0;
}
