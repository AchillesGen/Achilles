#include "nuchic/Integrators/BaseIntegrator.hh"

namespace nuchic {
namespace Integrator {

class AdaptiveIntegrator : public BaseIntegrator {
  public:
    AdaptiveIntegrator(bool cache = true) : BaseIntegrator(cache) {}
    AdaptiveIntegrator(const FunctionS &f, bool cache = true) : BaseIntegrator(f, cache) {}
    AdaptiveIntegrator(const FunctionV &f, bool cache = true) : BaseIntegrator(f, cache) {}
    AdaptiveIntegrator(const AdaptiveIntegrator &) = default;
    AdaptiveIntegrator(AdaptiveIntegrator &&) = default;
    AdaptiveIntegrator &operator=(const AdaptiveIntegrator &) = default;
    AdaptiveIntegrator &operator=(AdaptiveIntegrator &&) = default;

    ~AdaptiveIntegrator() override = default;

    virtual double Integrate(const double &, const double &, const double &, const double &) = 0;
    virtual std::vector<double> IntegrateVec(const double &, const double &, const double &,
                                             const double &) = 0;
};

} // namespace Integrator
} // namespace nuchic
