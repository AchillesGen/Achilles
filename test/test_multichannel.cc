#include "catch2/catch.hpp"
#include "nuchic/MultiChannel.hh"
#include "catch_utils.hh"

constexpr double s0 = -10.0;
constexpr double s1 = 10.0;

double test_func_exp(const std::vector<double> &x, double wgt) {
    double sms0 = x[0] - s0;
    double sms1 = x[0] - s1;
    return (2.0 * std::exp(-sms0 * sms0) + std::exp(-sms1 * sms1))
        / std::sqrt(std::acos(-1.0)) / 3.0 * wgt;
}

class DoubleMapper : public nuchic::Mapper<double> {
    public:
        DoubleMapper(size_t channel) : m_channel{std::move(channel)} {}
        void GeneratePoint(std::vector<double> &point, const std::vector<double> &rans) const override {
            double s = std::tan(std::acos(-1.0) * (rans[0] - 0.5));
            if(m_channel == 0) s += s0;
            else s += s1;
            point[0] = s;
        }
        double GenerateWeight(const std::vector<double> &point, std::vector<double> &rans) const override {
            const double sms0 = point[0] - s0;
            const double sms1 = point[0] - s1;
            if(m_channel == 1) {
                rans[0] = std::atan(sms1)/std::acos(-1.0) + 0.5;
                return 1.0 / (1.0 + sms1 * sms1) / std::acos(-1.0);
            }
            rans[0] = std::atan(sms0)/std::acos(-1.0) + 0.5;
            return 1.0 / (1.0 + sms0 * sms0) / std::acos(-1.0);
        }
    private:
        size_t m_channel;
};

TEST_CASE("YAML encoding / decoding Multichannel Summary", "[multichannel]") {
    nuchic::MultiChannelSummary summary;
    constexpr size_t nentries = 2;
    for(size_t i = 0; i < nentries; ++i) {
        auto vals = GENERATE(take(10, randomVector(100)));
        summary.results.emplace_back();
        for(const auto &val : vals) {
            summary.results.back() += val;
            summary.sum_results += val;
        }
    }
    auto weight = GENERATE(take(10, random(0.0, 1.0)));
    summary.best_weights.resize(2);
    summary.best_weights[0] = weight;
    summary.best_weights[1] = 1 - weight;

    YAML::Node node;
    node["Summary"] = summary;
    auto summary2 = node["Summary"].as<nuchic::MultiChannelSummary>();

    CHECK(summary.results.size() == summary2.results.size());
    for(size_t i = 0; i < summary.results.size(); ++i) 
        CHECK(summary.results[i] == summary2.results[i]);
    CHECK(summary.Result() == summary2.Result());

    CHECK(summary.best_weights == summary2.best_weights);
}

TEST_CASE("YAML encoding / decoding Multichannel Parameters", "[multichannel]") {
    nuchic::MultiChannelParams params{};

    YAML::Node node;
    node["Params"] = params;
    auto params2 = node["Params"].as<nuchic::MultiChannelParams>();

    CHECK(params.ncalls == params2.ncalls);
    CHECK(params.niterations == params2.niterations);
    CHECK(params.atol == params2.atol);
    CHECK(params.rtol == params2.rtol);
    CHECK(params.nrefine == params2.nrefine);
    CHECK(params.beta == params2.beta);
    CHECK(params.min_alpha == params2.min_alpha);
    CHECK(params.iteration == params2.iteration);
}

TEST_CASE("Multi-Channel Integration", "[multichannel]") {
    nuchic::Integrand<double> integrand(test_func_exp);
    for(size_t i = 0; i < 2; ++i) {
        nuchic::Channel<double> channel;
        nuchic::AdaptiveMap2 map(1, 50);
        channel.integrator = nuchic::Vegas2(map, nuchic::VegasParams{});
        channel.mapping = std::make_unique<DoubleMapper>(i);
        integrand.AddChannel(std::move(channel));
    }

    SECTION("Runs at least minimum required iterations") {
        static constexpr size_t nitn_min = 10;
        static constexpr double rtol = 2e-2, atol = 2e-2;

        nuchic::MultiChannel integrator(1, integrand.NChannels(),
                                        nuchic::MultiChannelParams{100000, nitn_min, rtol, atol});
        integrator.Optimize(integrand);
        auto results = integrator.Summary();

        CHECK(std::abs(results.sum_results.Mean() - 1.0) < 2*results.sum_results.Error());
        CHECK(results.results.size() >= nitn_min);
    }


    SECTION("Runs to desired precision") {
        static constexpr size_t nitn_min = 2;
        static constexpr double rtol = 1e-3, atol = 1e-3;
        nuchic::MultiChannel integrator(1, integrand.NChannels(),
                                        nuchic::MultiChannelParams{100000, nitn_min, rtol, atol});
        integrator.Optimize(integrand);
        auto results = integrator.Summary();

        CHECK(std::abs(results.sum_results.Mean() - 1.0) < 2*results.sum_results.Error());
        CHECK(results.sum_results.Error() < atol);
        CHECK(results.sum_results.Error()/results.sum_results.Mean() < rtol);
    }
}

TEST_CASE("YAML encoding / decoding Multichannel", "[multichannel]") {
    nuchic::Integrand<double> integrand(test_func_exp);
    for(size_t i = 0; i < 2; ++i) {
        nuchic::Channel<double> channel;
        nuchic::AdaptiveMap2 map(1, 50);
        channel.integrator = nuchic::Vegas2(map, nuchic::VegasParams{});
        channel.mapping = std::make_unique<DoubleMapper>(i);
        integrand.AddChannel(std::move(channel));
    }
    static constexpr size_t ncalls = 100000, nitn_min = 2;
    static constexpr double rtol = 1, atol = 1;
    nuchic::MultiChannel integrator(1, integrand.NChannels(),
                                    nuchic::MultiChannelParams{ncalls, nitn_min, rtol, atol});
    integrator.Optimize(integrand);
    auto results1 = integrator.Summary();

    YAML::Node node;
    node["Multichannel"] = integrator;

    auto integrator2 = node["Multichannel"].as<nuchic::MultiChannel>();
    auto results2 = integrator2.Summary();
    auto params = integrator2.Parameters();
   
    // Check integrator
    CHECK(integrator.Dimensions() == integrator2.Dimensions());
    CHECK(integrator.NChannels() == integrator2.NChannels());

    // Check parameters
    CHECK(params.ncalls == ncalls);
    CHECK(params.rtol == rtol);
    CHECK(params.atol == atol);
    CHECK(params.niterations == nitn_min);

    // Check previous result history
    CHECK(results1.sum_results.Mean() == results2.sum_results.Mean());
    CHECK(results1.sum_results.Error() == results2.sum_results.Error());
    CHECK(results1.best_weights == results2.best_weights);
}
