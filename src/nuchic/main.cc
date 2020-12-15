#include "nuchic/EventGen.hh"
#include "nuchic/Version.hh"
#include "nuchic/System.hh"
#include "nuchic/Logging.hh"

#include <dlfcn.h>

using namespace nuchic::SystemVariables;

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

int main() {

    Splash();
    CreateLogger(2, 5);

    const std::string path = "src/nuchic/fortran/";
    const std::string lib = libPrefix + "fortran_interface" + libSuffix;
    const std::string name = path + lib;
    void *handle = dlopen(name.c_str(), RTLD_NOW);
    if(!handle) {
        spdlog::warn("Cannot open HardScattering: {}", dlerror());
    }

    nuchic::EventGen generator("run.yml");
    generator.Initialize();
    generator.GenerateEvents();

    return 0;
}
