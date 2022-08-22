#ifndef MOM_SOLVER_HH
#define MOM_SOLVER_HH

#include <memory>
#include <utility>

namespace achilles {
    class FourVector;
    class Potential;

    FourVector SolveDelta(const FourVector&, const FourVector&, double, double, double, double);
    std::pair<double, double> FindMomentumRange(const FourVector&, double rangeExtend = 1.05);
    FourVector SolveDeltaWithPotential(const FourVector&, const Potential&,
                                       double, double, double, double, double, double);
}

#endif
