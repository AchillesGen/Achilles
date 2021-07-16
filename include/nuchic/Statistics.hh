#ifndef STATISTICS_HH
#define STATISTICS_HH

#include <algorithm>
#include <limits>

#include <iostream>
#include <cmath>

namespace nuchic {

using lim = std::numeric_limits<double>;

// Structure to hold moments
class StatsData {
    public:
        StatsData() = default;
        StatsData(const StatsData&) = default;
        StatsData(StatsData&&) = default;
        StatsData& operator=(const StatsData&) = default;
        StatsData& operator=(StatsData&&) = default;
        ~StatsData() = default;

        double Variance() const { return (sum2/n - Mean()*Mean()) / (n - 1); }

        StatsData& operator+=(double x) {
            n++;
            sum += x;
            sum2 += x*x;
            min = std::min(min, x);
            max = std::max(max, x);

            if(x != 0) n_finite++;

            return *this;
        }
        StatsData operator+(double x) {
            return {*this += x};
        }
        StatsData& operator+=(StatsData x) {
            n += x.n;
            n_finite += x.n_finite;
            sum += x.sum;
            sum2 += x.sum2;
            min = std::min(min, x.min);
            max = std::max(max, x.max);

            return *this;
        }
        StatsData operator+(const StatsData &x) {
            return {*this += x};
        }

        size_t Calls() const { return static_cast<size_t>(n); }
        size_t FiniteCalls() const { return static_cast<size_t>(n_finite); }
        double Mean() const { return sum/n; }
        double Min() const { return min; }
        double Max() const { return max; }
        double Error() const { return sqrt(Variance()); }

    private:
        double n{}, min{lim::max()}, max{lim::min()}, sum{}, sum2{}, n_finite{};
};

}

#endif
