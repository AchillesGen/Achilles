#include <algorithm>
#include <cmath>
#include <limits>

#include "nuchic/Integrators/GaussKronrod.hh"
#include "spdlog/spdlog.h"

using namespace nuchic::Integrators;

constexpr std::array<double, GaussKronrod::knots+1> GaussKronrod::KronrodWgts;
constexpr std::array<double, (GaussKronrod::knots+1)/2> GaussKronrod::GaussWgts;
constexpr std::array<double, GaussKronrod::knots> GaussKronrod::absc;

double GaussKronrod::Integrate(const double &a, const double &b, double &err, double&) {
    double halfLength = (b-a)*0.5; 
    double center = (a+b)*0.5;
    double fCenter = Function(center);

    double resultGauss = GaussWgts[GaussWgts.size()-1]*fCenter;
    double resultKronrod = KronrodWgts[KronrodWgts.size()-1]*fCenter;

    double absIntegral = std::abs(resultKronrod);

    for(size_t j = 1; j < KronrodWgts.size() - GaussWgts.size(); ++j) {
        size_t jj = 2*j - 1;
        double abscissa = halfLength * absc[jj];

        double f1 = Function(center - abscissa);
        double f2 = Function(center + abscissa);
        double fSum = f1 + f2;

        resultGauss += GaussWgts[j - 1] * fSum;
        resultKronrod += KronrodWgts[jj] * fSum;

        absIntegral += KronrodWgts[jj] * (std::abs(f1) + std::abs(f2));
    }

    for(size_t j = 0; j < GaussWgts.size(); ++j) {
        size_t jj = 2*j;
        double abscissa = halfLength * absc[jj];

        double f1 = Function(center - abscissa);
        double f2 = Function(center + abscissa);
        double fSum = f1 + f2;

        resultKronrod += KronrodWgts[jj] * fSum;

        absIntegral += KronrodWgts[jj] * (std::abs(f1) + std::abs(f2));
    }

    double resultMeanKronrod = resultKronrod * 0.5;
    double absDiffIntegral = KronrodWgts[KronrodWgts.size() - 1]
                           * (std::abs(fCenter - resultMeanKronrod));
    for(size_t j = 0; j < KronrodWgts.size(); ++j) {
        double abscissa = halfLength * absc[j];
        absDiffIntegral += std::abs(Function(center-abscissa) - resultMeanKronrod)
                         * KronrodWgts[j];
        absDiffIntegral += std::abs(Function(center+abscissa) - resultMeanKronrod)
                         * KronrodWgts[j];
    }

    absIntegral *= std::abs(halfLength);
    absDiffIntegral *= std::abs(halfLength);
    err = std::abs((resultKronrod - resultGauss)*halfLength);

    if(absDiffIntegral != 0 && err != 0)
        err = absDiffIntegral*std::min(1.0, pow(200*err/absDiffIntegral, 1.5));

    if(absIntegral > (std::numeric_limits<double>::min)() / (1E-12 * 50.))
        err = std::max(1E-12*50*absIntegral, err);

    spdlog::debug("GaussKronrod: I = {} +/- {}", resultKronrod*halfLength, err);
    return resultKronrod*halfLength;
}

std::vector<double> GaussKronrod::IntegrateVec(const double &a, const double &b,
                                               double &err, double&) {
    double halfLength = (b-a)*0.5; 
    double center = (a+b)*0.5;
    std::vector<double> fCenter = FunctionVec(center);

    auto resultGauss = fCenter;
    auto resultKronrod = fCenter;

    std::for_each(resultGauss.begin(), resultGauss.end(), [](double &x){
        x *= GaussWgts[GaussWgts.size()-1];
    });

    std::for_each(resultKronrod.begin(), resultKronrod.end(), [](double &x){
        x *= KronrodWgts[KronrodWgts.size()-1];
    });

    auto absIntegral = resultKronrod;
    std::for_each(absIntegral.begin(), absIntegral.end(), [](double &x){ x = std::abs(x); });

    for(size_t j = 1; j < KronrodWgts.size() - GaussWgts.size(); ++j) {
        size_t jj = 2*j - 1;
        double abscissa = halfLength * absc[jj];

        auto f1 = FunctionVec(center - abscissa);
        auto f2 = FunctionVec(center + abscissa);
        auto fSum = f1;
        std::transform(fSum.begin(), fSum.end(), f2.begin(), fSum.begin(), std::plus<double>());

        for(size_t i = 0; i < fSum.size(); ++i) {
            resultGauss[i] += GaussWgts[j - 1] * fSum[i];
            resultKronrod[i] += KronrodWgts[jj] * fSum[i];

            absIntegral[i] += KronrodWgts[jj] * (std::abs(f1[i]) + std::abs(f2[i]));
        }
    }

    for(size_t j = 0; j < GaussWgts.size(); ++j) {
        size_t jj = 2*j;
        double abscissa = halfLength * absc[jj];

        auto f1 = FunctionVec(center - abscissa);
        auto f2 = FunctionVec(center + abscissa);
        auto fSum = f1;
        std::transform(fSum.begin(), fSum.end(), f2.begin(), fSum.begin(), std::plus<double>());

        for(size_t i = 0; i < fSum.size(); ++i) {
            resultKronrod[i] += KronrodWgts[jj] * fSum[i];

            absIntegral[i] += KronrodWgts[jj] * (std::abs(f1[i]) + std::abs(f2[i]));
        }
    }

    auto resultMeanKronrod = resultKronrod;
    std::for_each(resultMeanKronrod.begin(), resultMeanKronrod.end(), [](double &x){ x *= 0.5; });
    std::vector<double> absDiffIntegral(resultMeanKronrod.size());
    std::transform(fCenter.begin(), fCenter.end(), resultMeanKronrod.begin(),
                   absDiffIntegral.begin(), [](const double &x1, const double &x2) {
                       return KronrodWgts[KronrodWgts.size() - 1] * std::abs(x1 - x2);
                   });

    for(size_t j = 0; j < KronrodWgts.size(); ++j) {
        double abscissa = halfLength * absc[j];
        auto f1 = FunctionVec(center-abscissa);
        auto f2 = FunctionVec(center+abscissa);
        for(size_t i = 0; i < f1.size(); ++i) {
            absDiffIntegral[i] += std::abs(f1[i] - resultMeanKronrod[i])
                             * KronrodWgts[j];
            absDiffIntegral[i] += std::abs(f2[i] - resultMeanKronrod[i])
                             * KronrodWgts[j];
        }
    }

    auto diff = resultKronrod;
    for(size_t i = 0; i < absIntegral.size(); ++i) {
        absIntegral[i] *= std::abs(halfLength);
        absDiffIntegral[i] *= std::abs(halfLength);
        diff[i] -= resultGauss[i];
        resultKronrod[i] *= halfLength;
    }
    err = std::abs(*std::max_element(diff.begin(), diff.end())*halfLength);

    auto absDiff = std::minmax_element(absDiffIntegral.begin(), absDiffIntegral.end());
    if(*absDiff.first != 0 && err != 0)
        err = *absDiff.second*std::min(1.0, pow(200*err/ *absDiff.first, 1.5));

    auto maxAbsIntegral = *std::max_element(absIntegral.begin(), absIntegral.end());
    if(maxAbsIntegral > (std::numeric_limits<double>::min)() / (1E-12 * 50.))
        err = std::max(1E-12*50*maxAbsIntegral, err);

    spdlog::debug("GaussKronrod: I = {} +/- {}", resultKronrod[0], err);
    return resultKronrod;
}
