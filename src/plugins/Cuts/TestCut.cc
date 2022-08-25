#include "Achilles/OneParticleCuts.hh"

namespace achilles {

ONE_PARTICLE_CUT(Test);

bool TestCut::MakeCut(const FourVector&) const {
    return true;
}

}
