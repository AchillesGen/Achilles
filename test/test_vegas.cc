#include "catch2/catch.hpp"

#include "nuchic/Vegas.hh"

double test_func(const std::vector<double> &x, double) {
    return 3.0/2.0*(x[0]*x[0] + x[1]*x[1]);
}

TEST_CASE("Integration", "[vegas]") {
    SECTION("Runs at least minimum required iterations") {
        static constexpr size_t nitn_min = 10;
        static constexpr double rtol = 1e-1, atol = 1e-1;
        nuchic::AdaptiveMap2 map(2, 100);
        nuchic::Vegas2 vegas(map, nuchic::VegasParams{10000, 5, rtol, atol, 1.5, nitn_min});

        vegas.Optimize(test_func);
        auto results = vegas.Summary();

        CHECK(std::abs(results.sum_results.Mean() - 1.0) < 2*results.sum_results.Error());
        CHECK(results.results.size() >= nitn_min);
    }

    SECTION("Runs to desired precision") {
        static constexpr size_t nitn_min = 2;
        static constexpr double rtol = 1e-3, atol = 1e-3;
        nuchic::AdaptiveMap2 map(2, 100);
        nuchic::Vegas2 vegas(map, nuchic::VegasParams{10000, 2, rtol, atol, 1.5, nitn_min});

        vegas.Optimize(test_func);
        auto results = vegas.Summary();

        CHECK(std::abs(results.sum_results.Mean() - 1.0) < 2*results.sum_results.Error());
        CHECK(results.sum_results.Error() < atol);
        CHECK(results.sum_results.Error()/results.sum_results.Mean() < rtol);
    }
}
