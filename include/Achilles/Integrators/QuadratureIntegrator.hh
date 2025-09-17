#include "Achilles/Integrators/BaseIntegrator.hh"

namespace achilles {
namespace Integrator {

class QuadratureIntegrator : public BaseIntegrator {
  public:
    QuadratureIntegrator(bool cache = true) : BaseIntegrator(cache) {}
    QuadratureIntegrator(const FunctionS &f, bool cache = true) : BaseIntegrator(f, cache) {}
    QuadratureIntegrator(const FunctionV &f, bool cache = true) : BaseIntegrator(f, cache) {}
    QuadratureIntegrator(const QuadratureIntegrator &) = default;
    QuadratureIntegrator(QuadratureIntegrator &&) = default;
    QuadratureIntegrator &operator=(const QuadratureIntegrator &) = default;
    QuadratureIntegrator &operator=(QuadratureIntegrator &&) = default;

    ~QuadratureIntegrator() override = default;

    virtual double Integrate(const double &, const double &, double &) = 0;
    virtual std::vector<double> IntegrateVec(const double &, const double &, double &) = 0;
};

} // namespace Integrator
} // namespace achilles
