#ifndef ADAPTIVE_INTEGRATOR_HH
#define ADAPTIVE_INTEGRATOR_HH

#include <memory>

#include "nuchic/Integrators/BaseIntegrator.hh"

namespace nuchic {

namespace Integrators {

class AdaptiveIntegrator : public BaseIntegrator {
    public:
        AdaptiveIntegrator(std::unique_ptr<BaseIntegrator>, const size_t&, bool cache=true);
        AdaptiveIntegrator(std::unique_ptr<BaseIntegrator>, const size_t&,
                           const FunctionD&, bool cache=true);
        AdaptiveIntegrator(std::unique_ptr<BaseIntegrator>, const size_t&,
                           const FunctionVD&, bool cache=true);

        void SetFunction(const FunctionD &f) { m_integrator -> SetFunction(f); Clear(); }
        void SetFunctionVec(const FunctionVD &fVec) { 
            m_integrator -> SetFunctionVec(fVec); 
            Clear() ;
        }

        double Integrate(const double&, const double&, double&, double&);
        std::vector<double> IntegrateVec(const double&, const double&, double&, double&);

    private:
        void Resize();
        void Clear();

        // Variables
        size_t m_maxSteps;
        std::unique_ptr<BaseIntegrator> m_integrator;
        std::vector<double> error, lowerList, upperList, resultList;
};

}

}

#endif
