#ifndef BASE_INTEGRATOR_HH
#define BASE_INTEGRATOR_HH

#include <functional>
#include <map>
#include <vector>

namespace nuchic {

namespace Integrators {

using FunctionD = std::function<double(const double&)>;
using FunctionVD = std::function<std::vector<double>(const double&)>;

class BaseIntegrator {
    public:
        BaseIntegrator(bool cache=true);
        BaseIntegrator(const FunctionD&, bool cache=true);
        BaseIntegrator(const FunctionVD&, bool cache=true);
        virtual ~BaseIntegrator() {}

        virtual void SetFunction(const FunctionD &f) {
            m_func = f;
            ClearCache();
        }
        virtual void SetFunctionVec(const FunctionVD &fVec) {
            m_funcVec = fVec;
            ClearCache();
        }
        void SetCache(const bool &cache) { m_cache = cache; }

        virtual double Integrate(const double&, const double&, double&, double&) = 0;
        virtual std::vector<double> IntegrateVec(const double&, const double&,
                                                 double&, double&) = 0;

    protected:
        double Function(const double&);
        std::vector<double> FunctionVec(const double&);
        void ClearCache();

    private:
        // Stored integrand
        FunctionD m_func;
        FunctionVD m_funcVec;

        // Function value cache
        bool m_cache;
        std::map<double, double> m_cacheFunc;
        std::map<double, std::vector<double>> m_cacheFuncVec;
};

}

}

#endif
