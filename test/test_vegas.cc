#include "catch2/catch.hpp"
#include "Achilles/Vegas.hh"
#include "catch_utils.hh"

double test_func(const std::vector<double> &x, double wgt) {
    return 3.0/2.0*(x[0]*x[0] + x[1]*x[1])*wgt;
}

double test_func2(const std::vector<double> &x, double wgt) {
    constexpr double s0 = -10.0;
    double s = std::tan(std::acos(-1.0) * (x[0] - 0.5)) + s0;

    double sms0 = s - s0;
    double wgt3 = 1.0 / (1.0 + sms0 * sms0) / std::acos(-1.0);
    // double wgt2 = std::acos(-1.0)/pow(cos(std::acos(-1.0)*(x[0] - 0.5)), 2);
    return (std::exp(-sms0 * sms0) )
          / std::sqrt(std::acos(-1.0)) * wgt / wgt3;

}

TEST_CASE("YAML encoding / decoding vegas summary", "[vegas]") {
    achilles::VegasSummary summary;
    constexpr size_t nentries = 3;
    for(size_t i = 0; i < nentries; ++i) {
        auto vals = GENERATE(take(10, randomVector(100)));
        summary.results.emplace_back();
        for(const auto &val : vals) {
            summary.results.back() += val;
            summary.sum_results += val;
        }
    }

    YAML::Node node;
    node["Summary"] = summary;
    auto summary2 = node["Summary"].as<achilles::VegasSummary>();

    CHECK(summary.results.size() == summary2.results.size());
    for(size_t i = 0; i < summary.results.size(); ++i)
        CHECK(summary.results[i] == summary2.results[i]);
    CHECK(summary.Result() == summary2.Result());
}

TEST_CASE("Vegas Integration", "[vegas]") {
    SECTION("Runs at least minimum required iterations") {
        static constexpr size_t nitn_min = 10;
        static constexpr double rtol = 2e-2, atol = 2e-2;
        achilles::AdaptiveMap map(2, 100);
        achilles::Vegas vegas(map, achilles::VegasParams{10000, 20, rtol, atol, 1.5, nitn_min});

        vegas.Optimize(test_func2);
        auto results = vegas.Summary();

        CHECK(std::abs(results.sum_results.Mean() - 1.0) < nsigma*results.sum_results.Error());
        CHECK(results.results.size() >= nitn_min);
    }

    SECTION("Runs to desired precision") {
        static constexpr size_t nitn_min = 2;
        static constexpr double rtol = 1e-3, atol = 1e-3;
        achilles::AdaptiveMap map(2, 100);
        achilles::Vegas vegas(map, achilles::VegasParams{10000, 2, rtol, atol, 1.5, nitn_min});

        vegas.Optimize(test_func);
        auto results = vegas.Summary();

        CHECK(std::abs(results.sum_results.Mean() - 1.0) < nsigma*results.sum_results.Error());
        CHECK(results.sum_results.Error() < atol);
        CHECK(results.sum_results.Error()/results.sum_results.Mean() < rtol);
    }
}

TEST_CASE("YAML encoding / decoding Vegas", "[vegas]") {
    static constexpr size_t nitn_min = 2;
    static constexpr double rtol = 1, atol = 1;
    achilles::AdaptiveMap map(2, 10);
    achilles::Vegas vegas(map, achilles::VegasParams{100, 2, rtol, atol, 1.5, nitn_min});
    vegas.Optimize(test_func);
    auto results1 = vegas.Summary();

    YAML::Node node;
    node["Vegas"] = vegas;

    auto vegas2 = node["Vegas"].as<achilles::Vegas>();
    auto results2 = vegas2.Summary();

    CHECK(results1.sum_results.Mean() == results2.sum_results.Mean());
    CHECK(results1.sum_results.Error() == results2.sum_results.Error());
}
