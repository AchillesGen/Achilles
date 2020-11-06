#include "nuchic/EventGen.hh"
#include "nuchic/System.hh"

#include <dlfcn.h>

using namespace nuchic::SystemVariables;

int main() {

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
