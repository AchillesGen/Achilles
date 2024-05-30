#ifndef STATISTICS_HH
#define STATISTICS_HH

#include <algorithm>
#include <iomanip>
#include <limits>

#include <cmath>
#include <iostream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

using lim = std::numeric_limits<double>;

// Class taken from: https://stackoverflow.com/q/3738349/9201027
class Percentile {
  public:
    Percentile(double percentile) : m_percentile(percentile) {}

    void Add(const double &x) {
        if(m_lower.empty() || x <= m_lower.front()) {
            m_lower.push_back(x);
            std::push_heap(m_lower.begin(), m_lower.end(), std::less<>());
        } else {
            m_upper.push_back(x);
            std::push_heap(m_upper.begin(), m_upper.end(), std::greater<>());
        }

        auto size = static_cast<size_t>(static_cast<double>(m_lower.size() + m_upper.size()) *
                                        m_percentile) +
                    1;
        if(m_lower.size() > size) {
            std::pop_heap(m_lower.begin(), m_lower.end(), std::less<>());
            m_upper.push_back(m_lower.back());
            std::push_heap(m_upper.begin(), m_upper.end(), std::greater<>());
            m_lower.pop_back();
        } else if(m_lower.size() < size) {
            std::pop_heap(m_upper.begin(), m_upper.end(), std::greater<>());
            m_lower.push_back(m_upper.back());
            std::push_heap(m_lower.begin(), m_lower.end(), std::less<>());
            m_upper.pop_back();
        }
    }

    double Get() const { return m_lower.front(); }

    void Clear() {
        m_lower.clear();
        m_upper.clear();
    }

    void SaveState(std::ostream &os) {
        const auto default_precision{os.precision()};
        os << std::setprecision(std::numeric_limits<double>::max_digits10 + 1);
        os << m_percentile << " ";
        os << m_lower.size() << " ";
        for(const auto &x : m_lower) os << x << " ";
        os << m_upper.size() << " ";
        for(const auto &x : m_upper) os << x << " ";
        os << std::setprecision(static_cast<int>(default_precision));
    }
    void LoadState(std::istream &is) {
        is >> m_percentile;
        size_t size;
        is >> size;
        m_lower.resize(size);
        for(auto &x : m_lower) is >> x;
        is >> size;
        m_upper.resize(size);
        for(auto &x : m_upper) is >> x;
    }

  private:
    double m_percentile;
    std::vector<double> m_lower, m_upper;
};

// Structure to hold moments
class StatsData {
  public:
    StatsData() = default;
    StatsData(const StatsData &) = default;
    StatsData(StatsData &&) = default;
    StatsData &operator=(const StatsData &) = default;
    StatsData &operator=(StatsData &&) = default;
    ~StatsData() = default;

    double Variance() const { return (sum2 / n - Mean() * Mean()) / (n - 1); }

    StatsData &operator+=(double x) {
        n++;
        sum += x;
        sum2 += x * x;
        min = std::min(min, x);
        max = std::max(max, x);

        if(x != 0) n_finite++;

        return *this;
    }
    StatsData operator+(double x) { return {*this += x}; }
    StatsData &operator+=(StatsData x) {
        n += x.n;
        n_finite += x.n_finite;
        sum += x.sum;
        sum2 += x.sum2;
        min = std::min(min, x.min);
        max = std::max(max, x.max);

        return *this;
    }
    StatsData operator+(const StatsData &x) { return {*this += x}; }

    size_t Calls() const { return static_cast<size_t>(n); }
    size_t FiniteCalls() const { return static_cast<size_t>(n_finite); }
    double Mean() const { return sum / n; }
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

    void SaveState(std::ostream &os) const {
        const auto default_precision{os.precision()};
        os << std::setprecision(std::numeric_limits<double>::max_digits10 + 1);
        os << n << " " << min << " " << max << " " << sum << " " << sum2 << " " << n_finite;
        os << std::setprecision(static_cast<int>(default_precision));
    }
    void LoadState(std::istream &is) { is >> n >> min >> max >> sum >> sum2 >> n_finite; }

    friend YAML::convert<achilles::StatsData>;

  private:
    double n{}, min{lim::max()}, max{lim::min()}, sum{}, sum2{}, n_finite{};
};

} // namespace achilles

namespace YAML {

template <> struct convert<achilles::StatsData> {
    static Node encode(const achilles::StatsData &rhs) {
        Node node;
        node = std::vector<double>{rhs.n, rhs.min, rhs.max, rhs.sum, rhs.sum2, rhs.n_finite};
        node.SetStyle(YAML::EmitterStyle::Flow);
        return node;
    }

    static bool decode(const Node &node, achilles::StatsData &rhs) {
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

} // namespace YAML

#endif
