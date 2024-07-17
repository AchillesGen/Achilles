#ifndef VEGAS_HH
#define VEGAS_HH

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Achilles/AdaptiveMap.hh"
#include "Achilles/Random.hh"
#include "Achilles/Statistics.hh"

#include "spdlog/spdlog.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

using lim = std::numeric_limits<double>;

class KBNSummation {
  public:
    KBNSummation() = default;
    inline double GetSum() noexcept { return sum + correction; }
    void AddTerm(double value) noexcept;
    inline void Reset() noexcept {
        sum = 0;
        correction = 0;
    }

  private:
    double sum{}, correction{};
};

template <typename T> using Func = std::function<double(const std::vector<T> &, const double &)>;

struct VegasParams {
    size_t ncalls{ncalls_default}, nrefine{nrefine_default};
    double rtol{rtol_default}, atol{atol_default}, alpha{alpha_default};
    size_t ninterations{nitn_default};

    static constexpr size_t nitn_default = 10, ncalls_default = 10000, nrefine_default = 5;
    static constexpr double alpha_default = 1.5, rtol_default = 1e-4, atol_default = 1e-4;
    static constexpr size_t nparams = 6;
};

struct VegasSummary {
    std::vector<StatsData> results;
    StatsData sum_results;

    StatsData Result() const { return sum_results; }
};

bool operator==(const VegasParams &lhs, const VegasParams &rhs);
bool operator==(const VegasSummary &lhs, const VegasSummary &rhs);

void SaveState(std::ostream &os, const VegasParams &params);
void LoadState(std::istream &is, VegasParams &params);
void SaveState(std::ostream &os, const VegasSummary &summary);
void LoadState(std::istream &is, VegasSummary &summary);

class Vegas {
  public:
    enum class Verbosity { silent, normal, verbose, very_verbose };

    Vegas() = default;
    Vegas(AdaptiveMap map, VegasParams _params)
        : grid{std::move(map)}, params{std::move(_params)} {}

    // Utilities
    void SetVerbosity(size_t v = 1) {
        if(v == 0)
            verbosity = Verbosity::silent;
        else if(v == 1)
            verbosity = Verbosity::normal;
        else if(v == 2)
            verbosity = Verbosity::verbose;
        else if(v == 3)
            verbosity = Verbosity::very_verbose;
        else
            throw std::runtime_error("Vegas: Invalid verbosity level");
    }
    AdaptiveMap Grid() const { return grid; }
    AdaptiveMap &Grid() { return grid; }

    bool operator==(const Vegas &rhs) const {
        return grid == rhs.grid && summary == rhs.summary && params == rhs.params;
    }
    void SaveState(std::ostream &os) const;
    void LoadState(std::istream &is);

    // Training the integratvegor
    void operator()(const Func<double> &);
    void Optimize(const Func<double> &);
    double GenerateWeight(const std::vector<double> &) const;
    void Adapt(const std::vector<double> &);
    void Refine();

    // Generating fixed number of events

    // Getting results
    VegasSummary Summary() const;

    // YAML interface
    friend YAML::convert<achilles::Vegas>;

  private:
    void PrintIteration() const;

    AdaptiveMap grid;
    VegasSummary summary;
    VegasParams params{};
    Verbosity verbosity{Verbosity::normal};
};

} // namespace achilles

namespace YAML {

template <> struct convert<achilles::VegasSummary> {
    static Node encode(const achilles::VegasSummary &rhs) {
        Node node;
        node["nentries"] = rhs.results.size();
        for(const auto &entry : rhs.results) node["entries"].push_back(entry);

        return node;
    }

    static bool decode(const Node &node, achilles::VegasSummary &rhs) {
        // Get the number of entries and ensure that is the number of entries
        // If the number of entries is zero, return then to prevent an error
        auto nentries = node["nentries"].as<size_t>();
        if(nentries == 0) return true;
        if(node["entries"].size() != nentries) return false;

        // Load the entries and keep track of the sum
        for(const auto &entry : node["entries"]) {
            rhs.results.push_back(entry.as<achilles::StatsData>());
            rhs.sum_results += rhs.results.back();
        }

        return true;
    }
};

template <> struct convert<achilles::Vegas> {
    static Node encode(const achilles::Vegas &rhs) {
        Node node;
        node["Grid"] = rhs.grid;
        node["Summary"] = rhs.summary;
        return node;
    }

    static bool decode(const Node &node, achilles::Vegas &rhs) {
        if(node.size() != 2) return false;

        rhs.grid = node["Grid"].as<achilles::AdaptiveMap>();
        rhs.summary = node["Summary"].as<achilles::VegasSummary>();
        return true;
    }
};

} // namespace YAML

#endif
