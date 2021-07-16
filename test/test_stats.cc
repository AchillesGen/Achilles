#include "catch2/catch.hpp"

#include "nuchic/Statistics.hh"

#include "catch_utils.hh"

TEST_CASE("Statistics class", "[vegas]") {
    SECTION("Adding individual points together") {
        nuchic::StatsData data;
        auto vals = GENERATE(take(100, randomVector(100)));
        double mean = 0, mean2 = 0, min = nuchic::lim::max(), max = nuchic::lim::min();
        for(const auto &val : vals) {
            data += val;
            mean += val;
            mean2 += val*val;
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
        CHECK(data.Variance() == Approx((mean2 - mean*mean)/static_cast<double>(vals.size()-1)));
    }

    SECTION("Adding multiple StatsData together") {
        nuchic::StatsData data1, data2;
        auto vals = GENERATE(take(100, randomVector(100)));
        double mean = 0, mean2 = 0, min = nuchic::lim::max(), max = nuchic::lim::min();
        bool odd = true;
        for(const auto &val : vals) {
            if(odd) data1 += val;
            else data2 += val;
            odd = !odd;
            mean += val;
            mean2 += val*val;
            if(val < min) min = val;
            if(val > max) max = val;
        }
        mean /= static_cast<double>(vals.size());
        mean2 /= static_cast<double>(vals.size());

        nuchic::StatsData data = data1 + data2;
        CHECK(data.Calls() == vals.size());
        CHECK(data.FiniteCalls() == vals.size());
        CHECK(data.Min() == min);
        CHECK(data.Max() == max);
        CHECK(data.Mean() == Approx(mean));
        CHECK(data.Variance() == Approx((mean2 - mean*mean)/static_cast<double>(vals.size()-1)));
    }
}
