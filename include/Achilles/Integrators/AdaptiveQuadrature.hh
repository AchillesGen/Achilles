#ifndef ADAPTIVEINTEGRATOR_HH
#define ADAPTIVEINTEGRATOR_HH

#include <memory>

#include "nuchic/Integrators/AdaptiveIntegrator.hh"

namespace nuchic {
namespace Integrator {

class QuadratureIntegrator;

class AdaptiveQuadrature : AdaptiveIntegrator {
  public:
    AdaptiveQuadrature(std::unique_ptr<QuadratureIntegrator>, const size_t &, bool cache = true);
    AdaptiveQuadrature(std::unique_ptr<QuadratureIntegrator>, const size_t &, const FunctionS &,
                       bool cache = true);
    AdaptiveQuadrature(std::unique_ptr<QuadratureIntegrator>, const size_t &, const FunctionV &,
                       bool cache = true);

    void SetFunction(const FunctionS &f) override;
    void SetFunctionVec(const FunctionV &fVec) override;

    double Integrate(const double &, const double &, const double &, const double &) override;
    std::vector<double> IntegrateVec(const double &, const double &, const double &,
                                     const double &) override;

  private:
    void Resize();
    void Clear();

    // Variables
    size_t m_maxSteps;
    std::unique_ptr<QuadratureIntegrator> m_integrator;
    std::vector<double> error, lowerList, upperList, resultList;
};

} // namespace Integrator
} // namespace nuchic

#endif
