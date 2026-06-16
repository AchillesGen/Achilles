#ifndef ACHILLES_FLOWINTEGRATOR_HH
#define ACHILLES_FLOWINTEGRATOR_HH

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "Achilles/FlowSampler.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Integrand.hh"
#include "Achilles/Integrator.hh"
#include "Achilles/Units.hh"

namespace achilles {

class Settings;

// Integrator backend that replaces the per-channel Vegas grids and channel-selection
// weights with a normalizing flow.
class FlowIntegrator : public Integrator {
  public:
    FlowIntegrator() = default;
    FlowIntegrator(size_t ndims, size_t nchannels, std::unique_ptr<FlowSampler> sampler,
                   UnitsEnum unit = UnitsEnum::nb, std::string cache_tag = "Flow");

    size_t Dimensions() const { return m_ndims; }
    size_t NChannels() const { return m_nchannels; }

    template <typename T> void operator()(Integrand<T> &);
    template <typename T> void Optimize(Integrand<T> &);
    template <typename T> std::vector<std::vector<T>> GeneratePoints(Integrand<T> &, size_t n);
    template <typename T>
    std::vector<double> GenerateWeights(Integrand<T> &, const std::vector<std::vector<T>> &points);

    // Integrator interface (runtime-polymorphic, FourVector only).
    void operator()(Integrand<FourVector> &) override;
    void Optimize(Integrand<FourVector> &) override;
    std::vector<std::vector<FourVector>> GeneratePoints(Integrand<FourVector> &, size_t n) override;
    std::vector<double>
    GenerateWeights(Integrand<FourVector> &,
                    const std::vector<std::vector<FourVector>> &points) override;
    void SetRtol(double rtol) override { m_rtol = rtol; }
    void SetNCalls(size_t ncalls) override { m_ncalls = ncalls; }

    bool HasResults() const override { return !m_results.empty(); }
    StatsData LastResult() const override { return m_results.back(); }
    bool NeedsOptimization(double rel_err) const override;
    bool RequiresTraining() const override { return m_requires_training; }
    std::string CacheTag() const override { return m_cache_tag; }
    void UpdateSummary(bool update) override { m_update_summary = update; }

    void SaveState(std::ostream &) const override;
    bool LoadState(std::istream &) override;

    // Tuning knobs (mirrors MultiChannelParams subset the interface needs).
    void SetNIterations(size_t n) { m_niterations = n; }
    void SetMaxIterations(size_t n) { m_max_iterations = n; }
    void SetWeightsPath(std::string path) { m_weights_path = std::move(path); }

  private:
    template <typename T>
    std::vector<double> MultiChannelWeights(Integrand<T> &func,
                                            const std::vector<std::vector<T>> &points);
    void PrintIteration() const;

    size_t m_ndims{};
    size_t m_nchannels{};
    std::unique_ptr<FlowSampler> m_sampler;
    UnitsEnum m_units{UnitsEnum::nb};

    double m_rtol{rtol_default};
    size_t m_ncalls{ncalls_default};
    size_t m_niterations{niterations_default};
    size_t m_max_iterations{max_iterations_default};
    bool m_requires_training{true};
    bool m_update_summary{true};
    std::vector<StatsData> m_results;
    std::string m_weights_path{"achilles_flow.weights"};
    std::string m_cache_tag{"Flow"};

    static constexpr double rtol_default = 1e-2;
    static constexpr size_t ncalls_default = 10000;
    static constexpr size_t niterations_default = 7;
    static constexpr size_t max_iterations_default = 1000;
};

std::unique_ptr<Integrator> MakeFlowIntegrator(const Settings &config, size_t ndims,
                                               size_t nchannels, UnitsEnum unit);

template <typename T>
std::vector<double>
achilles::FlowIntegrator::MultiChannelWeights(Integrand<T> &func,
                                              const std::vector<std::vector<T>> &points) {
    const size_t nch = func.NChannels();
    const size_t np = points.size();

    std::vector<size_t> channels;
    std::vector<std::vector<double>> rans;
    std::vector<double> map_density(np * nch);
    channels.reserve(np * nch);
    rans.reserve(np * nch);

    // For each point, inverse-map it under every channel: this yields the channel's
    // mapper density p_map_i(x) and the rans the flow density is evaluated on.
    for(size_t k = 0; k < np; ++k) {
        for(size_t i = 0; i < nch; ++i) {
            std::vector<double> r;
            map_density[k * nch + i] = func.GetChannel(i).mapping->GenerateWeight(points[k], r);
            channels.push_back(i);
            rans.push_back(std::move(r));
        }
    }

    const auto flow_density = m_sampler->Density(channels, rans);

    std::vector<double> weights(np);
    for(size_t k = 0; k < np; ++k) {
        double denom = 0;
        for(size_t i = 0; i < nch; ++i)
            denom += map_density[k * nch + i] * flow_density[k * nch + i];
        weights[k] = denom == 0 ? 0 : 1.0 / denom;
    }
    return weights;
}

template <typename T>
std::vector<std::vector<T>> achilles::FlowIntegrator::GeneratePoints(Integrand<T> &func, size_t n) {
    auto batch = m_sampler->Sample(n);
    std::vector<std::vector<T>> points(n, std::vector<T>(m_ndims));
    for(size_t i = 0; i < n; ++i)
        func.GetChannel(batch.channels[i]).mapping->GeneratePoint(points[i], batch.rans[i]);
    return points;
}

template <typename T>
std::vector<double>
achilles::FlowIntegrator::GenerateWeights(Integrand<T> &func,
                                          const std::vector<std::vector<T>> &points) {
    return MultiChannelWeights(func, points);
}

template <typename T> void achilles::FlowIntegrator::operator()(Integrand<T> &func) {
    // Draw the whole iteration in one flow forward pass.
    auto batch = m_sampler->Sample(m_ncalls);

    std::vector<std::vector<T>> points(m_ncalls, std::vector<T>(m_ndims));
    for(size_t i = 0; i < m_ncalls; ++i)
        func.GetChannel(batch.channels[i]).mapping->GeneratePoint(points[i], batch.rans[i]);

    // Batched multichannel phase-space weights (one Density call).
    const auto ps_wgts = MultiChannelWeights(func, points);

    StatsData results;
    std::vector<double> values(m_ncalls);
    for(size_t i = 0; i < m_ncalls; ++i) {
        const double val = ps_wgts[i] == 0 ? 0 : func(points[i], ps_wgts[i]);
        values[i] = val;
        results += val;
    }

    // One gradient step on the batch.
    m_sampler->Train(batch.channels, batch.rans, values);

    if(m_update_summary) m_results.push_back(results);
}

template <typename T> void achilles::FlowIntegrator::Optimize(Integrand<T> &func) {
    if(!m_requires_training) {
        // Pre-trained / inference-only flow (e.g. loaded from ONNX): the sampler will not
        // improve, so do not run the rtol-driven training loop. Just take m_niterations
        // passes to estimate the integral under the fixed distribution.
        for(size_t i = 0; i < std::max<size_t>(1, m_niterations); ++i) {
            (*this)(func);
            PrintIteration();
        }
        return;
    }

    double rel_err = std::numeric_limits<double>::max();
    while(NeedsOptimization(rel_err) && m_results.size() < m_max_iterations) {
        (*this)(func);
        const auto current = m_results.back();
        rel_err = current.Error() / std::abs(current.Mean());
        PrintIteration();
    }
}

} // namespace achilles

#endif
