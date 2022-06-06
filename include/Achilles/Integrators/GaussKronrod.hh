#ifndef GAUSSKRONROD_HH
#define GAUSSKRONROD_HH

#include "nuchic/Integrators/QuadratureIntegrator.hh"

namespace nuchic {
namespace Integrator {

class GaussKronrod : public QuadratureIntegrator {
    public:
        GaussKronrod(bool cache=true) : QuadratureIntegrator(cache) {}
        GaussKronrod(const FunctionS &func, bool cache=true) 
            : QuadratureIntegrator(func, cache) {}
        GaussKronrod(const FunctionV &func, bool cache)
            : QuadratureIntegrator(func, cache) {}

        double Integrate(const double&, const double&, double&) override;
        std::vector<double> IntegrateVec(const double&, const double&,
                                         double&) override;

    private:
        // Variables
        static constexpr size_t knots = 7;
        static const std::vector<double> KronrodWgts, GaussWgts, absc;
};

} // end Integrator namespace
} // end nuchic namespace

#endif
