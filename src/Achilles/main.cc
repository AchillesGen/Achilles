#include "Achilles/EventGen.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/Version.hh"
#include "Achilles/System.hh"
#include "Achilles/Logging.hh"
#ifdef ENABLE_BSM
#include "plugins/Sherpa/SherpaMEs.hh"
#include "plugins/Sherpa/Channels.hh"
#endif // ENABLE_BSM

#include "docopt.h"

#include <dlfcn.h>

using namespace achilles::SystemVariables;
using namespace achilles::PathVariables;

void Splash() {

    fmt::print(R"splash(
+=================================================================+
|                                                                 | 
|       d8b   db db    db  .o88b. db   db d888888b  .o88b.        |
|       888o  88 88    88 d8P  Y8 88   88   `88'   d8P  Y8        |
|       88V8o 88 88    88 8P      88ooo88    88    8P             |
|       88 V8o88 88    88 8b      88~~~88    88    8b             |
|       88  V888 88b  d88 Y8b  d8 88   88   .88.   Y8b  d8        |
|       VP   V8P ~Y8888P'  `Y88P' YP   YP Y888888P  `Y88P'        |
|                                                                 |
+-----------------------------------------------------------------+
|                                                                 |
|    Version: {:52}|
|    Authors: Joshua Isaacson, William Jay, Alessandro Lovato,    | 
|             Pedro A. Machado, Noemi Rocco                       | 
|                                                                 |
+=================================================================+
)splash", ACHILLES_VERSION);
}

static const std::string USAGE =
R"(
    Usage:
      achilles [<input>] [-v | -vv] [-s | --sherpa=<sherpa>...]
      achilles --display-cuts
      achilles --display-ps
      achilles --display-ff
      achilles (-h | --help)
      achilles --version

    Options:
      -v[v]                                 Increase verbosity level.
      -h --help                             Show this screen.
      --version                             Show version.
      -s <sherpa> --sherpa=<sherpa>         Define Sherpa option.
      --display-cuts                        Display the available cuts
      --display-ps                          Display the available phase spaces
      --display-ff                          Display the available form factors
)";

void GenerateEvents(const std::string &runcard, const std::vector<std::string> &shargs) {
    achilles::EventGen generator(runcard, shargs);
    generator.Initialize();
    generator.GenerateEvents();
}

int main(int argc, char *argv[]) {

    Splash();
    std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
                                                    { argv + 1, argv + argc },
                                                    true, // show help if requested
                                                    fmt::format("achilles {}", ACHILLES_VERSION)); //version string

    if(args["--display-cuts"].asBool()) {
        achilles::CutFactory<achilles::OneParticleCut>::DisplayCuts();
        achilles::CutFactory<achilles::TwoParticleCut>::DisplayCuts();
        return 0;
    }

    if(args["--display-ps"].asBool()) {
        achilles::PSFactory<achilles::HadronicBeamMapper, size_t>::DisplayPhaseSpaces();
        achilles::PSFactory<achilles::FinalStateMapper, std::vector<double>>::DisplayPhaseSpaces();
#ifdef ENABLE_BSM
        achilles::PSFactory<PHASIC::Channels, std::vector<double>>::DisplayPhaseSpaces();
#endif // ENABLE_BSM
        return 0;
    }

    if(args["--display-ff"].asBool()) {
        achilles::FormFactorFactory::Display();
        return 0;
    }

    std::string runcard = "run.yml";
    if(args["<input>"].isString()) runcard = args["<input>"].asString();
    
    auto verbosity = static_cast<int>(2 - args["-v"].asLong());
    CreateLogger(verbosity, 5);

    const std::string lib = libPrefix + "fortran_interface" + libSuffix;
    std::string name = installLibs + lib;
    void *handle = dlopen(name.c_str(), RTLD_NOW);
    if(!handle) {
        name = buildLibs + lib;
        handle = dlopen(name.c_str(), RTLD_NOW);
        if(!handle) {
            spdlog::warn("Cannot open HardScattering: {}", dlerror());
        }
    }
    std::vector<std::string> shargs;
    if (args["--sherpa"].isStringList()) shargs=args["--sherpa"].asStringList();

    GenerateEvents(runcard, shargs);

    // Close dynamic libraries
    if(handle) dlclose(handle);
    return 0;
}
