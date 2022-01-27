#ifndef VEGAS_HH
#define VEGAS_HH

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <future>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "nuchic/AdaptiveMap.hh"
#include "nuchic/Statistics.hh"
#include "nuchic/Random.hh"

#include "spdlog/spdlog.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace nuchic {

using lim = std::numeric_limits<double>;

class KBNSummation{
    public:
        KBNSummation() = default;
        inline double GetSum() noexcept { return sum + correction; }
        void AddTerm(double value) noexcept;
        inline void Reset() noexcept { sum = 0; correction = 0; }
    private:
        double sum{}, correction{};
};

using VegasPt = std::array<double,2>;
using Batch = std::vector<double>;
using BatchInt = std::vector<std::size_t>;
using Batch2D = std::vector<std::vector<double>>;
template<typename T>
using Func = std::function<double(const std::vector<T>&, const double&)>;

class VegasResult {
    public:
        VegasResult(bool weighted_ = true) : weighted(weighted_) {}
        void Update(double, double);
        VegasPt GetResult() const {return {mean, std::sqrt(var)};}
        bool Converged(const double&, const double&) const;
        std::string Summary() const;
        size_t Calls() const { return calls; }
        size_t FiniteCalls() const { return finiteCalls; }
        size_t NonZeroCalls() const { return nonZeroCalls; }
        double Max() const { return max; }

    private:
        size_t calls{}, nonZeroCalls{}, finiteCalls{};
        std::vector<double> sMean, sVar;
        double sum{}, sum2{}, sumVar{}, max{};
        double mean{}, meanDivVar{}, mean2DivVar{};
        double var{}, OneDivVar{};
        std::size_t n{};
        bool weighted;
        double chi2() const;
        double dof() const { return static_cast<double>(sMean.size() - 1); }
};

class Vegas {
    public:
        Vegas() = default;
        Vegas(nuchic::AdaptiveMap, const YAML::Node&);
        Vegas(const Vegas&) = default;
        Vegas(Vegas&&) = default;
        Vegas &operator=(const Vegas&) = default;
        Vegas &operator=(Vegas&&) = default;
        ~Vegas() = default;

        void SetDefaults();
        void Set(const YAML::Node&);
        void Clear() {result = VegasResult(adapt);}
        std::string Settings(const std::size_t& ngrid=0);
        void RandomBatch(std::size_t, std::size_t, Batch2D&, Batch2D&, Batch&, BatchInt&);
        VegasPt operator ()(const Func<double>&);
        
    private:
        static void SynchronizeMPI();

        void Setup();
//        void FillSigF(const std::function<double(std::vector<double>,double)>&,
//                      size_t cubeBase, size_t cubeBatch);
        Batch MakeBatch(const Func<double>&, Batch2D, Batch);
        Batch2D GenerateRandom(const std::size_t&, const std::size_t&);

        std::vector<double> sigF;
        Batch jac, fdv2;
        std::vector<std::size_t> nEvalCube;
        double sumSigF{};
        nuchic::AdaptiveMap map;
        std::size_t neval{}, nitn{}, maxinc{}, nCubeBatch{};
        std::size_t maxCube{}, maxEvalCube{}, minEvalCube{};
        std::size_t ndims{}, nstrats{}, nincrements{}, nCube{}, lastNEval{};
        double alpha{}, beta{}, rtol{}, atol{};
        bool adapt{}, adaptErrors{}, sync{};
        VegasResult result;

        // Default parameters
        static constexpr size_t nitn_default = 10;
        static constexpr size_t neval_default = 1000;
        static constexpr size_t maxinc_default = 1000;
        static constexpr size_t nCubeBatch_default = 1000;
        static constexpr size_t maxCube_default = 1e9;
        static constexpr size_t maxEvalCube_default = 1e7;
        static constexpr double alpha_default = 1.5;
        static constexpr double beta_default = 0.75;
        static constexpr bool adapt_default = true;
        static constexpr bool adaptErrors_default = false;
        static constexpr double rtol_default = 0;
        static constexpr double atol_default = 0;
        static constexpr bool sync_default = true;

// #ifdef USING_MPI
//         int GetMPIRank();
//         int mpiRank{};
// #endif
};

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

class Vegas2 {
    public:
        enum class Verbosity {
            silent,
            normal,
            verbose,
            very_verbose
        };

        Vegas2() = default;
        Vegas2(AdaptiveMap2 map, VegasParams _params) : grid{std::move(map)}, params{std::move(_params)} {}

        // Utilities
        void SetVerbosity(size_t v = 1) {
            if(v == 0) verbosity = Verbosity::silent; 
            else if(v == 1) verbosity = Verbosity::normal;
            else if(v == 2) verbosity = Verbosity::verbose;
            else if(v == 3) verbosity = Verbosity::very_verbose;
            else throw std::runtime_error("Vegas: Invalid verbosity level");
        }
        AdaptiveMap2 Grid() const { return grid; }
        AdaptiveMap2 &Grid() { return grid; }
        // bool Serialize(std::ostream &out) const {
        //     
        // }

        // Training the integratvegor
        void operator()(const Func<double>&);
        void Optimize(const Func<double>&);
        double GenerateWeight(const std::vector<double>&) const;
        void Adapt(const std::vector<double>&);
        void Refine();

        // Generating fixed number of events

        // Getting results
        VegasSummary Summary() const; 

        // YAML interface
        friend YAML::convert<nuchic::Vegas2>;

    private:
        void PrintIteration() const;

        AdaptiveMap2 grid;
        VegasSummary summary;
        VegasParams params{};
        Verbosity verbosity{Verbosity::normal};
};

}

namespace YAML {

template<>
struct convert<nuchic::VegasSummary> {
    static Node encode(const nuchic::VegasSummary &rhs) {
        Node node;
        node["nentries"] = rhs.results.size();
        for(const auto &entry : rhs.results)
            node["entries"].push_back(entry);

        return node;
    }

    static bool decode(const Node &node, nuchic::VegasSummary &rhs) {
        // Get the number of entries and ensure that is the number of entries
        // If the number of entries is zero, return then to prevent an error
        auto nentries = node["nentries"].as<size_t>();
        if(nentries == 0) return true;
        if(node["entries"].size() != nentries) return false;

        // Load the entries and keep track of the sum
        for(const auto &entry : node["entries"]) {
            rhs.results.push_back(entry.as<nuchic::StatsData>());
            rhs.sum_results += rhs.results.back();
        }

        return true;
    }
};

template<>
struct convert<nuchic::Vegas2> {
    static Node encode(const nuchic::Vegas2 &rhs) {
        Node node;
        node["Grid"] = rhs.grid;
        node["Summary"] = rhs.summary;
        return node;
    }

    static bool decode(const Node &node, nuchic::Vegas2 &rhs) {
        if(node.size() != 2) return false;

        rhs.grid = node["Grid"].as<nuchic::AdaptiveMap2>();
        rhs.summary = node["Summary"].as<nuchic::VegasSummary>();
        return true;
    }
};

}

#endif
