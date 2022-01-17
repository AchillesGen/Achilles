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
        inline double GetSum() noexcept {return sum + correction;};
        void AddTerm(double value) noexcept;
        inline void Reset() noexcept {sum = 0; correction = 0;} ;
    private:
        double sum{}, correction{};
};

using VegasPt = std::array<double,2>;
using Batch = std::vector<double>;
using BatchInt = std::vector<std::size_t>;
using Batch2D = std::vector<std::vector<double>>;
using Func = std::function<double(const std::vector<double>&, const double&)>;

class VegasResult {
    public:
        VegasResult(bool weighted_ = true) : weighted(weighted_) {}
        void Update(double, double);
        VegasPt GetResult() const {return {mean, std::sqrt(var)};}
        bool Converged(const double&, const double&) const;
        std::string Summary() const;

    private:
        std::vector<double> sMean, sVar;
        double sum{}, sum2{}, sumVar{};
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
        VegasPt operator ()(const Func&);
        
    private:
        static void SynchronizeMPI();

        void Setup();
//        void FillSigF(const std::function<double(std::vector<double>,double)>&,
//                      size_t cubeBase, size_t cubeBatch);
        Batch MakeBatch(const Func&, Batch2D, Batch);
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

}

#endif
