#ifndef DOUBLEEXPONENTIAL_HH
#define DOUBLEEXPONENTIAL_HH

#include <array>
#include <cmath>

#include "Achilles/Integrators/AdaptiveIntegrator.hh"
#include "Achilles/Utilities.hh"

namespace achilles {
namespace Integrator {

// Setup abscissa and weights for Double Exponential integration
constexpr size_t _phases = 6;
constexpr size_t _size = 6 * 1 << _phases;

struct DEPoint {
    double abscissa;
    double weights;
};

inline double t2(const size_t &k) {
    return exp(static_cast<double>(k) / ipow(2, _phases));
}

inline double u1(const size_t &k) {
    return M_PI_4 * (t2(k) + 1.0 / t2(k));
}

inline double t3(const size_t &k) {
    return exp(M_PI_4 * (t2(k) - 1.0 / t2(k)));
}

inline double t4(const size_t &k) {
    return (t3(k) + 1.0 / t3(k)) / 2;
}

inline double GetAbcs(const size_t &k) {
    return 1.0 / (t3(k) * t4(k));
}

inline double GetWeight(const size_t &k) {
    return u1(k) / (t4(k) * t4(k));
}

class DoubleExponential : public AdaptiveIntegrator {
  public:
    DoubleExponential(bool cache = false) : AdaptiveIntegrator(cache) {}
    DoubleExponential(const FunctionS &func, bool cache = false)
        : AdaptiveIntegrator(func, cache) {}
    DoubleExponential(const FunctionV &func, bool cache = false)
        : AdaptiveIntegrator(func, cache) {}

    double Integrate(const double &, const double &, const double &, const double &) override;
    std::vector<double> IntegrateVec(const double &, const double &, const double &,
                                     const double &) override;

  private:
    // Table for Double Exponential integration
    static std::array<DEPoint, _size> table;
    static bool initialized;

    DEPoint GeneratePoint(size_t curr) { return {GetAbcs(curr), GetWeight(curr)}; }
};

} // namespace Integrator
} // namespace achilles

#endif
