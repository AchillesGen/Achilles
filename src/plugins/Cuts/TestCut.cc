#include "Achilles/OneParticleCuts.hh"
#include "Achilles/Version.hh"

EXPORT_ACHILLES_VERSION()

namespace achilles {

ONE_PARTICLE_CUT(Test);

bool TestCut::MakeCut(const FourVector&) const {
    return true;
}

void __attribute__((destructor)) CleanUp() {
    spdlog::debug("Destructor called!");
    CutFactory<OneParticleCut>::Deregister("Test");
}

}
