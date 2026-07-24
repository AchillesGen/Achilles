#include "catch2/catch.hpp"

#include "Achilles/Statistics.hh"

#include "catch_utils.hh"

TEST_CASE("Statistics class", "[vegas]") {
    SECTION("Adding individual points together") {
        achilles::StatsData data;
        auto vals = GENERATE(take(100, randomVector(100)));
        double mean = 0, mean2 = 0, min = achilles::lim::max(), max = achilles::lim::min();
        for(const auto &val : vals) {
            data += val;
            mean += val;
            mean2 += val * val;
            if(val < min) min = val;
            if(val > max) max = val;
        }
        mean /= static_cast<double>(vals.size());
        mean2 /= static_cast<double>(vals.size());

        CHECK(data.Calls() == vals.size());
        CHECK(data.FiniteCalls() == vals.size());
        CHECK(data.Min() == min);
        CHECK(data.Max() == max);
        CHECK(data.Mean() == Approx(mean));
        CHECK(data.Variance() ==
              Approx((mean2 - mean * mean) / static_cast<double>(vals.size() - 1)));
    }

    SECTION("Adding multiple StatsData together") {
        achilles::StatsData data1, data2;
        auto vals = GENERATE(take(100, randomVector(100)));
        double mean = 0, mean2 = 0, min = achilles::lim::max(), max = achilles::lim::min();
        bool odd = true;
        for(const auto &val : vals) {
            if(odd)
                data1 += val;
            else
                data2 += val;
            odd = !odd;
            mean += val;
            mean2 += val * val;
            if(val < min) min = val;
            if(val > max) max = val;
        }
        mean /= static_cast<double>(vals.size());
        mean2 /= static_cast<double>(vals.size());

        achilles::StatsData data = data1 + data2;
        CHECK(data.Calls() == vals.size());
        CHECK(data.FiniteCalls() == vals.size());
        CHECK(data.Min() == min);
        CHECK(data.Max() == max);
        CHECK(data.Mean() == Approx(mean));
        CHECK(data.Variance() ==
              Approx((mean2 - mean * mean) / static_cast<double>(vals.size() - 1)));
    }
}

TEST_CASE("YAML encoding / decoding StatsData", "[vegas]") {
    achilles::StatsData data1, data2;
    auto vals = GENERATE(take(100, randomVector(100)));
    for(const auto &val : vals) { data1 += val; }

    YAML::Node node;
    node["Stats"] = data1;
    data2 = node["Stats"].as<achilles::StatsData>();

    CHECK(data1.Calls() == data2.Calls());
    CHECK(data1.Min() == data2.Min());
    CHECK(data1.Max() == data2.Max());
    CHECK(data1.Mean() == data2.Mean());
    CHECK(data1.Error() == data2.Error());
    CHECK(data1.FiniteCalls() == data2.FiniteCalls());
}

TEST_CASE("Rejected (zero-weight) trials count toward the mean", "[vegas]") {
    // The geometry-driver sigma estimator is sum(w) / N_trials with one entry
    // per generator trial; rejected trials enter as zeros and must dilute the
    // mean, while FiniteCalls tracks only the emitted events.
    achilles::StatsData data;
    data += 4.0;
    data += 0.0;
    data += 0.0;
    data += 0.0;

    CHECK(data.Calls() == 4);
    CHECK(data.FiniteCalls() == 1);
    CHECK(data.Mean() == Approx(1.0));
}

TEST_CASE("Variance is finite for degenerate samples", "[vegas]") {
    // A single entry has no defined spread; it must not divide by zero.
    achilles::StatsData single;
    single += 2.5;
    CHECK(single.Variance() == 0.0);
    CHECK(single.Error() == 0.0);

    // Identical entries: rounding can push sum2/n - mean^2 slightly negative,
    // which must not turn into a NaN error via sqrt.
    achilles::StatsData equal;
    for(size_t i = 0; i < 1000; ++i) equal += 0.1;
    CHECK(equal.Variance() >= 0.0);
    CHECK(std::isfinite(equal.Error()));
}
