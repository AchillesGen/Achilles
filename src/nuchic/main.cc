#include "nuchic/EventGen.hh"
#include "nuchic/Version.hh"
#include "nuchic/System.hh"
#include "nuchic/Logging.hh"
#include "plugins/SherpaMEs.hh"

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
      nuchic (-h | --help)
      nuchic --version

    Options:
      -v[v]                                 Increase verbosity level.
      -h --help                             Show this screen.
      --version                             Show version.
      -s <sherpa> --sherpa=<sherpa>         Define Sherpa option.
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
    for (auto arg: shargs) std::cout<<arg<<std::endl;
    sherpa.Initialize(shargs);

    GenerateEvents(runcard,&sherpa);

    // Close dynamic libraries
    dlclose(handle);
    return 0;
}
