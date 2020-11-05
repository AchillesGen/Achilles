#include "nuchic/EventGen.hh"

#include <dlfcn.h>

int main() {

    void *handle = dlopen("src/nuchic/fortran/libfortran_interface.so", RTLD_NOW);
    if(!handle) {
        spdlog::warn("Cannot open HardScattering: {}", dlerror());
    }

    nuchic::EventGen generator("run.yml");
    generator.Initialize();
    generator.GenerateEvents();

    return 0;
}
