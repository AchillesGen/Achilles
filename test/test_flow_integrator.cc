#include "Achilles/FlowIntegrator.hh"
#include "Achilles/Random.hh"
#include "catch2/catch.hpp"
#include "catch_utils.hh"

#include <sstream>

// Same analytic target as the MultiChannel test: two shifted Gaussians whose integral
// over the real line is 1, importance-sampled with two Cauchy channel maps.
constexpr double s0 = -10.0;
constexpr double s1 = 10.0;

static double test_func_exp(const std::vector<double> &x, double wgt) {
    double sms0 = x[0] - s0;
    double sms1 = x[0] - s1;
    return (2.0 * std::exp(-sms0 * sms0) + std::exp(-sms1 * sms1)) / std::sqrt(std::acos(-1.0)) /
           3.0 * wgt;
}

class DoubleMapper : public achilles::Mapper<double> {
  public:
    DoubleMapper(size_t channel) : m_channel{channel} {}
    void GeneratePoint(std::vector<double> &point, const std::vector<double> &rans) override {
        double s = std::tan(std::acos(-1.0) * (rans[0] - 0.5));
        s += (m_channel == 0) ? s0 : s1;
        point[0] = s;
    }
    double GenerateWeight(const std::vector<double> &point, std::vector<double> &rans) override {
        const double sms0 = point[0] - s0;
        const double sms1 = point[0] - s1;
        rans.resize(1);
        if(m_channel == 1) {
            rans[0] = std::atan(sms1) / std::acos(-1.0) + 0.5;
            return 1.0 / (1.0 + sms1 * sms1) / std::acos(-1.0);
        }
        rans[0] = std::atan(sms0) / std::acos(-1.0) + 0.5;
        return 1.0 / (1.0 + sms0 * sms0) / std::acos(-1.0);
    }
    size_t NDims() const override { return 1; }

  private:
    size_t m_channel;
};

// A non-learning sampler that draws channels and unit-cube points uniformly. With it the
// FlowIntegrator reduces to plain multichannel Monte Carlo (uniform channel weights, no
// learned grid), so it must still converge to the known integral. It also lets us assert
// the integrator drives the sampler each iteration.
class MockFlowSampler : public achilles::FlowSampler {
  public:
    MockFlowSampler(size_t ndims, size_t nchannels)
        : m_ndims{ndims}, m_density{1.0 / static_cast<double>(nchannels)},
          m_weights(nchannels, 1.0) {}

    achilles::FlowSampleBatch Sample(size_t n) override {
        ++m_sample_calls;
        achilles::FlowSampleBatch batch;
        batch.channels.resize(n);
        batch.rans.assign(n, std::vector<double>(m_ndims));
        batch.density.assign(n, m_density);
        for(size_t i = 0; i < n; ++i) {
            achilles::Random::Instance().Generate(batch.rans[i]);
            batch.channels[i] = achilles::Random::Instance().SelectIndex(m_weights);
        }
        return batch;
    }
    std::vector<double> Density(const std::vector<size_t> &channels,
                                const std::vector<std::vector<double>> &) override {
        return std::vector<double>(channels.size(), m_density);
    }
    double Train(const std::vector<size_t> &, const std::vector<std::vector<double>> &,
                 const std::vector<double> &) override {
        ++m_train_calls;
        return 0.0;
    }
    void Save(const std::string &) const override {}
    bool Load(const std::string &) override { return false; }

    size_t SampleCalls() const { return m_sample_calls; }
    size_t TrainCalls() const { return m_train_calls; }

  private:
    size_t m_ndims;
    double m_density;
    std::vector<double> m_weights;
    size_t m_sample_calls{};
    size_t m_train_calls{};
};

static achilles::Integrand<double> MakeIntegrand() {
    achilles::Integrand<double> integrand(test_func_exp);
    for(size_t i = 0; i < 2; ++i) {
        achilles::Channel<double> channel;
        channel.mapping = std::make_unique<DoubleMapper>(i);
        achilles::AdaptiveMap map(channel.mapping->NDims(), 2);
        channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
        integrand.AddChannel(std::move(channel));
    }
    return integrand;
}

TEST_CASE("FlowIntegrator converges with a uniform sampler", "[flow]") {
    auto integrand = MakeIntegrand();
    auto sampler = std::make_unique<MockFlowSampler>(1, integrand.NChannels());
    auto *sampler_ptr = sampler.get();

    constexpr size_t niterations = 5;
    achilles::FlowIntegrator flow(1, integrand.NChannels(), std::move(sampler));
    flow.SetNCalls(20000);
    flow.SetNIterations(niterations);
    flow.SetRtol(1.0); // satisfiable; the non-learning sampler keeps per-iteration error flat

    flow.Optimize(integrand);

    // The integrator drove the sampler each iteration: one Sample()/Train() per iteration
    // (plus the batched weight queries via Density, not counted here).
    CHECK(sampler_ptr->TrainCalls() == niterations);
    CHECK(sampler_ptr->SampleCalls() == niterations);
    REQUIRE(flow.HasResults());

    const auto result = flow.LastResult();
    CHECK(std::abs(result.Mean() - 1.0) < nsigma * result.Error());
}

TEST_CASE("FlowIntegrator GeneratePoints/GenerateWeights are batched and finite", "[flow]") {
    auto integrand = MakeIntegrand();
    auto sampler = std::make_unique<MockFlowSampler>(1, integrand.NChannels());
    achilles::FlowIntegrator flow(1, integrand.NChannels(), std::move(sampler));

    constexpr size_t n = 256;
    auto points = flow.GeneratePoints(integrand, n);
    auto weights = flow.GenerateWeights(integrand, points);

    REQUIRE(points.size() == n);
    REQUIRE(weights.size() == n);
    for(size_t i = 0; i < n; ++i) {
        CHECK(points[i].size() == 1);
        CHECK(std::isfinite(weights[i]));
        CHECK(weights[i] > 0);
    }
}

TEST_CASE("FlowIntegrator cache rejects a foreign backend tag", "[flow]") {
    // A MultiChannel-tagged stream must not be consumed by the flow LoadState; it should
    // report a mismatch (false) and stay un-optimized.
    auto sampler = std::make_unique<MockFlowSampler>(1, 2);
    achilles::FlowIntegrator flow(1, 2, std::move(sampler));

    std::stringstream ss;
    ss << static_cast<int>(achilles::IntegratorType::MultiChannel) << " 1 2 ";
    CHECK_FALSE(flow.LoadState(ss));
    CHECK_FALSE(flow.HasResults());
}
