#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdlib>

#include "nuchic/Integrators/AdaptiveIntegrator.hh"
#include "spdlog/spdlog.h"

using namespace nuchic::Integrators;

AdaptiveIntegrator::AdaptiveIntegrator(std::unique_ptr<BaseIntegrator> integrator,
                                       const size_t &maxSteps, bool cache) 
    : BaseIntegrator(cache),  m_maxSteps(maxSteps), m_integrator(std::move(integrator)) {
    
    m_integrator -> SetCache(cache);    
    Resize();
}

AdaptiveIntegrator::AdaptiveIntegrator(std::unique_ptr<BaseIntegrator> integrator, 
                                       const size_t &maxSteps, const FunctionD &func, bool cache)
    : BaseIntegrator(func, cache), m_maxSteps(maxSteps), m_integrator(std::move(integrator)) {

    m_integrator -> SetFunction(func);
    m_integrator -> SetCache(cache);    
    Resize();
}

AdaptiveIntegrator::AdaptiveIntegrator(std::unique_ptr<BaseIntegrator> integrator, 
                                       const size_t &maxSteps, const FunctionVD &func, bool cache)
    : BaseIntegrator(func, cache), m_maxSteps(maxSteps), m_integrator(std::move(integrator)) {

    m_integrator -> SetFunctionVec(func);
    m_integrator -> SetCache(cache);    
    Resize();
}

void AdaptiveIntegrator::Resize() {
    error.resize(m_maxSteps, 0);
    lowerList.resize(m_maxSteps, 0);
    upperList.resize(m_maxSteps, 0);
    resultList.resize(m_maxSteps, 0);
}

void AdaptiveIntegrator::Clear() {
    std::fill(error.begin(), error.end(), 0);
    std::fill(lowerList.begin(), lowerList.end(), 0);
    std::fill(upperList.begin(), upperList.end(), 0);
    std::fill(resultList.begin(), resultList.end(), 0);
}

std::vector<double> AdaptiveIntegrator::IntegrateVec(const double&, const double&,
                                                     double&, double&) {
    throw std::runtime_error("Not implemented!");
}

double AdaptiveIntegrator::Integrate(const double &a, const double &b,
                                     double &rerr, double &aerr) {
    double err, dummy;
    double result = m_integrator -> Integrate(a, b, err, dummy);
    error[0] = err;
    lowerList[0] = a;
    upperList[0] = b;
    resultList[0] = result;

    double errorBound = std::max(aerr, rerr*std::abs(result));

    spdlog::trace("Iteration 1:");
    spdlog::trace("  - lower = {}, upper = {}, value = {}, error = {}",
                  lowerList[0], upperList[0], result, error[0]);

    if(err < errorBound) return result;

    size_t maxError = 0;
    double area = result;
    double errorSum = err;
    error[maxError] = err;

    size_t nSubintervals;

    for(nSubintervals = 2; nSubintervals <= m_maxSteps; ++nSubintervals) {
        double lower1 = lowerList[maxError];
        double upper2 = upperList[maxError];

        double upper1 = (lower1+upper2)*0.5;
        double lower2 = upper1;

        double error1, error2;

        double area1 = m_integrator -> Integrate(lower1, upper1, error1, dummy);
        double area2 = m_integrator -> Integrate(lower2, upper2, error2, dummy);

        spdlog::trace("Iteration {}:", nSubintervals);
        spdlog::trace("  - lower1 = {}, upper1 = {}, value1 = {}, error1 = {}",
                      lower1, upper1, area1, error1);
        spdlog::trace("  - lower2 = {}, upper2 = {}, value2 = {}, error2 = {}",
                      lower2, upper2, area2, error2);

        double area12 = area1 + area2;
        double error12 = error1 + error2;
        errorSum += error12 - error[maxError];
        area += area12 - resultList[maxError];

        if(error2 > error1) {
            lowerList[maxError] = lower2;
            lowerList[nSubintervals-1] = lower1;
            upperList[maxError] = upper2;
            upperList[nSubintervals-1] = upper1;
            resultList[maxError] = area2;
            resultList[nSubintervals-1] = area1;
            error[maxError] = error2;
            error[nSubintervals-1] = error1;
        } else {
            lowerList[nSubintervals-1] = lower2;
            upperList[maxError] = upper1;
            upperList[nSubintervals-1] = upper2;
            upperList[maxError] = upper1;
            resultList[nSubintervals-1] = area2;
            resultList[maxError] = area1;
            error[nSubintervals-1] = error2;
            error[maxError] = error1;
        }

        maxError = static_cast<size_t>(std::distance(error.begin(),
            std::max_element(error.begin(), error.end())));

        errorBound = std::max(aerr, rerr*std::abs(area));

        if(errorSum <= errorBound) break;
    }

    spdlog::trace("Integration used {} subintervals", nSubintervals);
    return std::accumulate(resultList.begin(),
                           resultList.begin()+static_cast<int>(nSubintervals), 0.0);
}
