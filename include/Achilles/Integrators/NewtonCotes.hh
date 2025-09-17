#ifndef NEWTONCOTES_HH
#define NEWTONCOTES_HH

#include "Achilles/Integrators/QuadratureIntegrator.hh"

namespace achilles {
namespace Integrator {

class NewtonCotes : public QuadratureIntegrator {
  public:
    NewtonCotes(const size_t &order, bool closed = true);
    NewtonCotes(const size_t &order, bool closed = true, bool cache = true);
    NewtonCotes(const size_t &order, const FunctionS &func, bool closed = true, bool cache = true);
    NewtonCotes(const size_t &order, const FunctionV &func, bool closed = true, bool cache = true);

    double Integrate(const double &, const double &, double &) override;
    std::vector<double> IntegrateVec(const double &, const double &, double &) override;

  private:
    void SetupClosed();
    void SetupOpen();

    size_t m_order;
    std::function<double(const double &, const double &)> baseFunc, errFunc;
    FunctionV baseFuncV, errFuncV;
};

} // namespace Integrator
} // namespace achilles

#endif
