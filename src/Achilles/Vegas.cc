#include <algorithm>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <stdexcept>
#include <utility>

#include "Achilles/AdaptiveMap.hh"
#include "Achilles/Random.hh"
#include "Achilles/Vegas.hh"

#ifdef USING_MPI
// #include "Achilles/MPI.hh"
#endif

void achilles::Vegas::operator()(const Func<double> &func) {
    std::vector<double> rans(grid.Dims());
    std::vector<double> train_data(grid.Dims() * grid.Bins());

    StatsData results;

    for(size_t i = 0; i < params.ncalls; ++i) {
        Random::Instance().Generate(rans);

        double wgt = grid(rans);
        double val = func(rans, wgt);
        double val2 = val * val;

        results += val;

        for(size_t j = 0; j < grid.Dims(); ++j) {
            train_data[j * grid.Bins() + grid.FindBin(j, rans[j])] += val2;
        }
    }

    grid.Adapt(params.alpha, train_data);
    summary.results.push_back(results);
    summary.sum_results += results;
}

void achilles::Vegas::Optimize(const Func<double> &func) {
    double abs_err = lim::max(), rel_err = lim::max();
    size_t irefine = 0;
    while((abs_err > params.atol && rel_err > params.rtol) ||
          summary.results.size() < params.ninterations) {
        (*this)(func);
        StatsData current = summary.Result();
        abs_err = current.Error();
        rel_err = abs_err / std::abs(current.Mean());

        PrintIteration();
        if(++irefine == params.nrefine) {
            Refine();
            irefine = 0;
        }
    }
}

double achilles::Vegas::GenerateWeight(const std::vector<double> &rans) const {
    return grid.GenerateWeight(rans);
}

void achilles::Vegas::Adapt(const std::vector<double> &train_data) {
    grid.Adapt(params.alpha, train_data);
}

void achilles::Vegas::Refine() {
    grid.Split();
    params.ncalls *= 2;
}

achilles::VegasSummary achilles::Vegas::Summary() const {
    std::cout << "Final integral = "
              << fmt::format("{:^8.5e} +/- {:^8.5e} ({:^8.5e} %)", summary.Result().Mean(),
                             summary.Result().Error(),
                             summary.Result().Error() / summary.Result().Mean() * 100)
              << std::endl;
    return summary;
}

void achilles::Vegas::PrintIteration() const {
    std::cout << fmt::format("{:3d}   {:^8.5e} +/- {:^8.5e}    {:^8.5e} +/- {:^8.5e}",
                             summary.results.size(), summary.results.back().Mean(),
                             summary.results.back().Error(), summary.Result().Mean(),
                             summary.Result().Error())
              << std::endl;
}

// Below are functions for preforming the modified Kahan Summation, to ensure
// that the number of threads does not effect the result for a fixed random seed
void achilles::KBNSummation::AddTerm(double value) noexcept {
    /// Function to add a value to the sum that is being calculated and keep the
    /// correction term
    double t = sum + value;
    if(std::abs(sum) >= std::abs(value)) {
        correction += ((sum - t) + value);
    } else {
        correction += ((value - t) + sum);
    }
    sum = t;
}
