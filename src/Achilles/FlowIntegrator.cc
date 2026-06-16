#include "Achilles/FlowIntegrator.hh"

#include <iomanip>
#include <limits>

#include "spdlog/spdlog.h"

using achilles::FlowIntegrator;

FlowIntegrator::FlowIntegrator(size_t ndims, size_t nchannels, std::unique_ptr<FlowSampler> sampler,
                               UnitsEnum unit, std::string cache_tag)
    : m_ndims{ndims}, m_nchannels{nchannels}, m_sampler{std::move(sampler)}, m_units{unit},
      m_cache_tag{std::move(cache_tag)} {
    if(m_sampler) m_requires_training = m_sampler->RequiresTraining();
}

bool FlowIntegrator::NeedsOptimization(double rel_err) const {
    return (rel_err > m_rtol) || m_results.size() < m_niterations;
}

void FlowIntegrator::PrintIteration() const {
    const auto &current = m_results.back();
    const double scale = UnitScale(m_units);
    const std::string unit = ToString(m_units);
    spdlog::info("{:3d}   {:^8.5e} +/- {:^8.5e} {} ( {:^4.2g} %)", m_results.size(),
                 current.Mean() * scale, current.Error() * scale, unit,
                 current.Error() / current.Mean() * 100);
}

// Integrator interface: forward to the templated cores instantiated on FourVector.
void FlowIntegrator::operator()(Integrand<FourVector> &func) {
    this->operator()<FourVector>(func);
}

void FlowIntegrator::Optimize(Integrand<FourVector> &func) {
    this->Optimize<FourVector>(func);
}

std::vector<std::vector<achilles::FourVector>>
FlowIntegrator::GeneratePoints(Integrand<FourVector> &func, size_t n) {
    return this->GeneratePoints<FourVector>(func, n);
}

std::vector<double>
FlowIntegrator::GenerateWeights(Integrand<FourVector> &func,
                                const std::vector<std::vector<FourVector>> &points) {
    return this->GenerateWeights<FourVector>(func, points);
}

void FlowIntegrator::SaveState(std::ostream &os) const {
    const auto default_precision{os.precision()};
    os << std::setprecision(std::numeric_limits<double>::max_digits10 + 1);
    os << static_cast<int>(IntegratorType::Flow) << " ";
    os << m_ndims << " " << m_nchannels << " ";
    os << m_rtol << " " << m_ncalls << " " << m_niterations << " ";
    os << m_results.size() << " ";
    for(const auto &result : m_results) {
        result.SaveState(os);
        os << " ";
    }
    os << std::setprecision(static_cast<int>(default_precision));

    // The learned flow parameters live outside the text stream.
    if(m_sampler) m_sampler->Save(m_weights_path);
}

bool FlowIntegrator::LoadState(std::istream &is) {
    // A missing cache, or one written by a different integrator backend, is not an error:
    // leave this integrator empty so it reports NeedsOptimization() and re-trains. Do not
    // consume the rest of the stream in that case.
    int tag;
    if(!(is >> tag)) return false;
    if(tag != static_cast<int>(IntegratorType::Flow)) return false;

    is >> m_ndims >> m_nchannels;
    is >> m_rtol >> m_ncalls >> m_niterations;

    size_t nresults;
    is >> nresults;
    m_results.clear();
    m_results.resize(nresults);
    for(auto &result : m_results) result.LoadState(is);

    if(m_sampler && !m_sampler->Load(m_weights_path)) {
        // Summary restored but the trained weights could not be: drop the summary so the
        // flow is retrained rather than generating from an untrained network.
        spdlog::warn("FlowIntegrator: could not load flow weights from {}; will retrain.",
                     m_weights_path);
        m_results.clear();
        return false;
    }

    return true;
}
