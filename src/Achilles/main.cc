#include "Achilles/EventGen.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Logging.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Particle.hh"
#include "Achilles/System.hh"
#include "Achilles/Version.hh"
#include "Achilles/fortran/FNuclearModel.hh"
#include "git.h"
#include "plugins/Manager/PluginManager.hh"
#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/Channels.hh"
#include "plugins/Sherpa/SherpaInterface.hh"
#endif

#include "docopt.h"

#include <dlfcn.h>
#include <filesystem>

using namespace achilles::SystemVariables;
using namespace achilles::PathVariables;
namespace fs = std::filesystem;

void Splash() {
    fmt::print(R"splash(
+=====================================================================+
|                                                                     | 
|    .d8b.   .o88b. db   db d888888b db      db      d88888b .d8888.  |
|   d8' `8b d8P  Y8 88   88   `88'   88      88      88'     88'  YP  |
|   88ooo88 8P      88ooo88    88    88      88      88ooooo `8bo.    |
|   88~~~88 8b      88~~~88    88    88      88      88~~~~~  `Y8b.   |
|   88   88 Y8b  d8 88   88   .88.   88booo. 88booo. 88.     db  8D   |
|   YP   YP  `Y88P' YP   YP Y888888P Y88888P Y88888P Y88888P `8888Y'  |
|                                                                     | 
+---------------------------------------------------------------------+
|                                                                     |
|                 .d88888888888888888888888888b.                      |
|            .d8888888888888888888888888888888888b.                   |
|         .d88888888888888888888888888888888888888888b.               |
|      88888888888888888888888888888888888888888888888888b.           |
|       `Y8888888888888888888888888888888888888888888888888b.         |
|            `Y88888888888888888888888888888888888888888888888        |
|              `Y88888888888888888888888888888888888888888888888      |
|                 `Y888888888888888888888888888P  Y88888888888888     |
|                    `Y88  8888  8888  88Y'        Y8888888888888     |
|                   .d888888888888888888888b.       Y888888888888.    |
|                .d88888888888888888888888888b.      8888888888888    |
|              .d8888888888888888888888888888888b   88888888888888    |
|            .d8888888888888888888888888888888888   8888888888888D    |
|            d88888888888888888888888888888888888   8888888888888     |
|           d888888888888888888888888888888888888   Y888888888888     |
|          d888888P'   d888888888888888888888888P       88888888      |
|          8888P     .d888888888888888888888888P         888888       |
|          888P    .d8888888888888888888888888P          8888Y        |
|          888b  .d888888888888888888888888888           8888P        |
|          Y88  d888888888888888888888888888888           888         |
|          `8' d8888888888888888888888888888888            88         |
|              8888888888888888888888888888888P             8         |
|              8888888888888888P     Y8888888P                        |
|             d8888888888888888       Y88888P                         |
|             8888888888888888P        Y888P                          |
|            d8888888888888P            Y8P                           |
|            888888888P                  8                            |
|           d888P                                                     |
|                                                                     |
+---------------------------------------------------------------------+
|                                                                     |
|    Version: {:56}|
|    Authors: Joshua Isaacson, William Jay, Alessandro Lovato,        | 
|        Pedro A. Machado, Luke Pickering, Noemi Rocco,               | 
|        Noah Steinberg                                               |
|                                                                     |
|    Undergraduate Student Contributions:                             |
|        Diego Lopez Gutierrez, Sherry Wang, Russell Farnsworth       |
|                                                                     |
+=====================================================================+
)splash",
               ACHILLES_VERSION);
}

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
      achilles [<input>] [-v | -vv] [-l | -ll] [-s | --sherpa=<sherpa>...]
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
      -h --help                             Show this screen.
      --version                             Show version.
      -s <sherpa> --sherpa=<sherpa>         Define Sherpa option.
      --display-cuts                        Display the available cuts
      --display-ps                          Display the available phase spaces
      --display-ff                          Display the available form factors
      --display-int-models                  Display the available cascade interaction models
      --display-nuc-models                  Display the available nuclear interaction models
)";

void GenerateEvents(const std::string &runcard, const std::vector<std::string> &shargs) {
    achilles::EventGen generator(runcard, shargs);
    generator.Initialize();
    generator.GenerateEvents();
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

    auto verbosity = static_cast<int>(2 - args["-v"].asLong());
    auto log_verbosity = std::min(verbosity, static_cast<int>(2 - args["-l"].asLong()));
    CreateLogger(verbosity, log_verbosity, 1);
    GitInformation();
    achilles::Plugin::Manager plugin_manager;

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
    if(args["<input>"].isString())
        runcard = args["<input>"].asString();
    else {
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

    try {
        GenerateEvents(runcard, shargs);
    } catch(const std::runtime_error &error) { spdlog::error(error.what()); }

    return 0;
}
