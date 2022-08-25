#include "Achilles/OneParticleCuts.hh"
#include "Achilles/Version.hh"


#undef ACHILLES_VERSION_MINOR
#define ACHILLES_VERSION_MINOR 1

EXPORT_ACHILLES_VERSION()

namespace achilles {

ONE_PARTICLE_CUT(Test);

bool TestCut::MakeCut(const FourVector&) const {
    return true;
}

void __attribute__((destructor)) CleanUp() {
    spdlog::info("Destructor called!");
    CutFactory<OneParticleCut>::Deregister("Test");
}

}
