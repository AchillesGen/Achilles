#include "nuchic/Integrators/BaseIntegrator.hh"

using namespace nuchic::Integrator;

BaseIntegrator::BaseIntegrator(bool cache) : m_cache(cache) {
    m_func = nullptr;
    m_funcVec = nullptr;
}

BaseIntegrator::BaseIntegrator(FunctionS func, bool cache)
    : m_func(std::move(func)), m_cache(cache) {
    m_funcVec = nullptr;
}

BaseIntegrator::BaseIntegrator(FunctionV funcVec, bool cache)
    : m_funcVec(std::move(funcVec)), m_cache(cache) {
    m_func = nullptr;
}

double BaseIntegrator::Function(const double &x) {
    if(m_cache) {
        auto it = m_cacheFunc.find(x);
        double fx;
        if(it == m_cacheFunc.end()) {
            fx = m_func(x);
            m_cacheFunc[x] = fx;
        } else {
            fx = it->second;
        }
        return fx;
    }

    return m_func(x);
}

std::vector<double> BaseIntegrator::FunctionVec(const double &x) {
    if(m_cache) {
        auto it = m_cacheFuncVec.find(x);
        std::vector<double> fx;
        if(it == m_cacheFuncVec.end()) {
            fx = m_funcVec(x);
            m_cacheFuncVec[x] = fx;
        } else {
            fx = it->second;
        }
        return fx;
    }

    return m_funcVec(x);
}

void BaseIntegrator::ClearCache() {
    m_cacheFunc.clear();
    m_cacheFuncVec.clear();
}
