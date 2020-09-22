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

namespace nuchic {

using lim = std::numeric_limits<double>;

class KBNSummation{
    public:
        KBNSummation(): sum{0}, correction{0} { };
        inline double GetSum() noexcept {return sum + correction;};
        void AddTerm(double value) noexcept;
        inline void Reset() noexcept {sum = 0; correction = 0;} ;
    private:
        double sum, correction;
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
        double mean{}, mean2{}, meanDivVar{}, mean2DivVar{};
        double var{}, var2{}, OneDivVar{};
        std::size_t n{};
        bool weighted;
        double chi2() const;
        double dof() const { return static_cast<double>(sMean.size() - 1); }
};

class Vegas {
    public:
        Vegas(nuchic::AdaptiveMap, const std::map<std::string, std::string>&);
        Vegas(const Vegas&) = default;
        Vegas(Vegas&&) = default;
        Vegas &operator=(const Vegas&) = default;
        Vegas &operator=(Vegas&&) = default;
        ~Vegas() { delete rand; }

        void Set(const std::map<std::string,std::string>&);
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
        randutils::mt19937_rng* rand;

// #ifdef USING_MPI
//         int GetMPIRank();
//         int mpiRank{};
// #endif
};

}

#endif
