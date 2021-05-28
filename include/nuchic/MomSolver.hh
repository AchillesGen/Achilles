#ifndef MOM_SOLVER_HH
#define MOM_SOLVER_HH

namespace nuchic {
    class FourVector;

    FourVector SolveDelta(const FourVector&, const FourVector&, double, double, double, double);
    FourVector SolveDeltaWithPotential(const FourVector&, const FourVector&,
                                       double, double, double, double, double);
}

#endif
