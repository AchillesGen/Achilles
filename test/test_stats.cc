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
