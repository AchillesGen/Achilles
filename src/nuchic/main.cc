#include "nuchic/EventGen.hh"
#include "nuchic/FinalStateMapper.hh"
#include "nuchic/HadronicMapper.hh"
#include "nuchic/Version.hh"
#include "nuchic/System.hh"
#include "nuchic/Logging.hh"
#include "plugins/Sherpa/SherpaMEs.hh"
#include "plugins/Sherpa/Channels.hh"

#include "docopt.h"

#include <dlfcn.h>

using namespace nuchic::SystemVariables;
using namespace nuchic::PathVariables;

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
)splash", NUCHIC_VERSION);
}

static const std::string USAGE =
R"(
    Usage:
      nuchic [<input>] [-v | -vv] [-s | --sherpa=<sherpa>...]
      nuchic --display-cuts
      nuchic --display-ps
      nuchic (-h | --help)
      nuchic --version

    Options:
      -v[v]                                 Increase verbosity level.
      -h --help                             Show this screen.
      --version                             Show version.
      -s <sherpa> --sherpa=<sherpa>         Define Sherpa option.
      --display-cuts                        Display the available cuts
      --display-ps                          Display the available phase spaces
)";

void GenerateEvents(const std::string &runcard,nuchic::SherpaMEs *const sherpa) {
    nuchic::EventGen generator(runcard,sherpa);
    generator.Initialize();
    generator.GenerateEvents();
}

int main(int argc, char *argv[]) {

    Splash();
    std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
                                                    { argv + 1, argv + argc },
                                                    true, // show help if requested
                                                    fmt::format("Nuchic {}", NUCHIC_VERSION)); //version string

    if(args["--display-cuts"].asBool()) {
        nuchic::CutFactory<nuchic::OneParticleCut>::DisplayCuts();
        nuchic::CutFactory<nuchic::TwoParticleCut>::DisplayCuts();
        return 0;
    }

    if(args["--display-ps"].asBool()) {
        nuchic::PSFactory<nuchic::HadronicBeamMapper, size_t>::DisplayPhaseSpaces();
        nuchic::PSFactory<nuchic::FinalStateMapper, std::vector<double>>::DisplayPhaseSpaces();
        nuchic::PSFactory<PHASIC::Channels, std::vector<double>>::DisplayPhaseSpaces();
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
    nuchic::SherpaMEs sherpa;
    std::vector<std::string> shargs;
    if (args["--sherpa"].isStringList()) shargs=args["--sherpa"].asStringList();
    sherpa.Initialize(shargs);

    GenerateEvents(runcard,&sherpa);

    // Close dynamic libraries
    dlclose(handle);
    return 0;
}
