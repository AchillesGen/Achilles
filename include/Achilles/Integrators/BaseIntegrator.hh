#ifndef BASE_INTEGRATOR_HH
#define BASE_INTEGRATOR_HH

#include <functional>
#include <map>
#include <vector>

namespace nuchic {
namespace Integrator {

using FunctionS = std::function<double(const double &)>;
using FunctionV = std::function<std::vector<double>(const double &)>;

class BaseIntegrator {
  public:
    BaseIntegrator(bool cache = true);
    BaseIntegrator(FunctionS, bool cache = true);
    BaseIntegrator(FunctionV, bool cache = true);
    BaseIntegrator(const BaseIntegrator &) = default;
    BaseIntegrator(BaseIntegrator &&) = default;
    BaseIntegrator &operator=(const BaseIntegrator &) = default;
    BaseIntegrator &operator=(BaseIntegrator &&) = default;

    virtual ~BaseIntegrator() = default;

    virtual void SetFunction(const FunctionS &f) {
        m_func = f;
        ClearCache();
    }
    virtual void SetFunctionVec(const FunctionV &fVec) {
        m_funcVec = fVec;
        ClearCache();
    }
    void SetCache(const bool &cache) { m_cache = cache; }

  protected:
    double Function(const double &);
    std::vector<double> FunctionVec(const double &);
    void ClearCache();

  private:
    // Stored integrand
    FunctionS m_func;
    FunctionV m_funcVec;

    // Function value cache
    bool m_cache;
    std::map<double, double> m_cacheFunc;
    std::map<double, std::vector<double>> m_cacheFuncVec;
};

} // namespace Integrator
} // namespace nuchic

#endif
