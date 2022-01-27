#ifndef STATISTICS_HH
#define STATISTICS_HH

#include <algorithm>
#include <limits>

#include <iostream>
#include <cmath>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

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

        bool operator==(const StatsData &other) const {
            static constexpr double tol = 1e-6;
            bool equal = n == other.n && n_finite == other.n_finite;
            equal = equal && (std::abs(min - other.min) < tol);
            equal = equal && (std::abs(max - other.max) < tol);
            equal = equal && (std::abs(sum - other.sum) < tol);
            equal = equal && (std::abs(sum2 - other.sum2) < tol);
            return equal;
        }
        bool operator!=(const StatsData &other) const { return !(*this == other); }

        friend YAML::convert<nuchic::StatsData>;

    private:
        double n{}, min{lim::max()}, max{lim::min()}, sum{}, sum2{}, n_finite{};
};

}

namespace YAML {

template<>
struct convert<nuchic::StatsData> {
    static Node encode(const nuchic::StatsData &rhs) {
        Node node;
        node = std::vector<double>{rhs.n, rhs.min, rhs.max, rhs.sum, rhs.sum2, rhs.n_finite};
        node.SetStyle(YAML::EmitterStyle::Flow);
        return node;
    }

    static bool decode(const Node &node, nuchic::StatsData &rhs) {
        // Ensure the node has 6 entries
        if(node.size() != 6 || !node.IsSequence()) return false;

        // Load the entries
        rhs.n = node[0].as<double>();
        rhs.min = node[1].as<double>();
        rhs.max = node[2].as<double>();
        rhs.sum = node[3].as<double>();
        rhs.sum2 = node[4].as<double>();
        rhs.n_finite = node[5].as<double>();

        return true;
    }
};

}

#endif
