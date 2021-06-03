#ifndef MOM_SOLVER_HH
#define MOM_SOLVER_HH

#include <utility>

namespace nuchic {
    class FourVector;

    FourVector SolveDelta(const FourVector&, const FourVector&, double, double, double, double);
    std::pair<double, double> FindMomentumRange(double, double, double, double, double);
    FourVector SolveDeltaWithPotential(const FourVector&,
                                       double, double, double, double, double, double);
    double Potential(double p, double rho);
}

#endif
