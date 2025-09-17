#include "Achilles/Integrators/DoubleExponential.hh"
#include "spdlog/spdlog.h"

using namespace achilles::Integrator;

bool DoubleExponential::initialized = false;
std::array<DEPoint, _size> DoubleExponential::table = {};

double DoubleExponential::Integrate(const double &x1, const double &x2, const double &ceps1,
                                    const double &ceps2) {
    if(!initialized) {
        size_t idx = 0;
        for(auto &row : table) row = GeneratePoint(idx++);
        initialized = true;
    }

    static constexpr int izx = 5;
    double ax = (x2 - x1) / 2;
    double bx = (x2 + x1) / 2;
    double err = 0.0, h = 1.0;
    double t1, t2, tw1, tw2;

    static constexpr size_t phases = _phases + 1;
    static constexpr std::array<size_t, phases> ip = {1, 2, 4, 8, 16, 32, 64};

    size_t evals = 0;
    double sum = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    for(size_t k = 1; k <= _phases; ++k) {
        h /= 2;
        s3 = s2;
        s2 = s1;
        double fm1 = 0.0;
        double fm2 = 0.0;
        size_t k1 = ip[_phases - k];
        size_t k2 = ip[_phases - k + 1];

        size_t iz1 = 0, iz2 = 0;

        // Evaluate function at level k in x, avoiding unnecessary computation
        for(size_t i = 0; i < _size; i += k1) {
            if(i % k2 != 0 || k == 1) {
                double xt1 = 1.0 - table[i].abscissa;
                double xx1 = bx - ax * xt1;
                double xx2 = bx + ax * xt1;
                bool log1 = xx1 > x1 && iz1 < izx;
                bool log2 = xx2 < x2 && iz2 < izx;

                if(!log1 && !log2) break;

                if(log1) {
                    t1 = Function(xx1);
                    evals++;
                    tw1 = t1 * table[i].weights;
                    if(std::abs(tw1) < ceps1)
                        iz1++;
                    else
                        iz1 = 0;
                } else {
                    t1 = 0.0;
                    tw1 = 0.0;
                }

                if(i > 0 && log2) {
                    t2 = Function(xx2);
                    evals++;
                    tw2 = t2 * table[i].weights;
                    if(std::abs(tw2) < ceps1)
                        iz2++;
                    else
                        iz2 = 0;
                } else {
                    t2 = 0.0;
                    tw2 = 0.0;
                }

                sum += tw1 + tw2;
                tw1 = std::abs(tw1);
                tw2 = std::abs(tw2);
                t1 = std::abs(t1);
                t2 = std::abs(t2);

                fm1 = std::max({fm1, tw1, tw2});
                fm2 = std::max({fm2, t1, t2});
            }
        }

        s1 = ax * h * sum;
        double eps1 = fm1 * ceps1;
        double eps2 = fm2 * ceps2;
        double d1 = log10(std::abs(s1 - s2));
        double d2 = log10(std::abs(s1 - s3));
        double d3 = log10(eps1) - 1;
        double d4 = log10(eps2) - 1;

        if(k <= 2)
            err = 1.0;
        else if(d1 == -16)
            err = 0.0;
        else {
            double val1 = d1 * d1 / d2;
            double val2 = 2 * d1;
            double max = std::max({val1, val2, d3, d4});
            err = pow(10, static_cast<int>(std::min(0.0, max)));
        }

        spdlog::trace("Iteration {} (value, error): {}, {}", k, s1, err);
        spdlog::trace("Number of evaluations = {}", evals);

        if(k > 3 && err < eps1) {
            spdlog::trace("Result (value, error): {}, {}", s1, err);
            return s1;
        }
        if(k >= 3 && err < eps2) {
            spdlog::trace("Result (value, error): {}, {}", s1, err);
            return s1;
        }
    }

    if(s1 != 0)
        spdlog::warn("DoubleExponential: Failed to converge (value, error): {}, {}", s1, err);
    return s1;
}

std::vector<double> DoubleExponential::IntegrateVec(const double &x1, const double &x2,
                                                    const double &ceps1, const double &ceps2) {
    if(!initialized) {
        size_t idx = 0;
        for(auto &row : table) row = GeneratePoint(idx++);
    }

    static constexpr int izx = 5;
    double ax = (x2 - x1) / 2;
    double bx = (x2 + x1) / 2;
    double err = 0.0, h = 1.0;
    static size_t size = FunctionVec(x1).size();

    static constexpr size_t phases = _phases + 1;
    static constexpr std::array<size_t, phases> ip = {1, 2, 4, 8, 16, 32, 64};

    size_t evals = 0;
    std::vector<double> sum(size, 0.0), s1(size, 0.0), s2(size, 0.0), s3(size, 0.0);
    std::vector<double> t1(size), t2(size), tw1(size), tw2(size);
    for(size_t k = 1; k <= _phases; ++k) {
        h /= 2;
        std::vector<double> fm1(size, 0.0);
        std::vector<double> fm2(size, 0.0);
        s3 = s2;
        s2 = s1;
        size_t k1 = ip[_phases - k];
        size_t k2 = ip[_phases - k + 1];

        size_t iz1 = 0, iz2 = 0;

        // Evaluate function at level k in x, avoiding unnecessary computation
        for(size_t i = 0; i < _size; i += k1) {
            if(i % k2 != 0 || k == 1) {
                double xt1 = 1.0 - table[i].abscissa;
                double xx1 = bx - ax * xt1;
                double xx2 = bx + ax * xt1;
                bool log1 = xx1 > x1 && iz1 < izx;
                bool log2 = xx2 < x2 && iz2 < izx;

                if(!log1 && !log2) break;

                if(log1) {
                    t1 = FunctionVec(xx1);
                    evals++;
                    std::transform(t1.begin(), t1.end(), tw1.begin(),
                                   [&](const double x) { return x * table[i].weights; });
                    if(std::abs(*std::max_element(tw1.begin(), tw1.end())) < ceps1)
                        iz1++;
                    else
                        iz1 = 0;
                } else {
                    std::fill(t1.begin(), t1.end(), 0.0);
                    std::fill(tw1.begin(), tw1.end(), 0.0);
                }

                if(i > 0 && log2) {
                    t2 = FunctionVec(xx2);
                    evals++;
                    std::transform(t2.begin(), t2.end(), tw2.begin(),
                                   [&](const double x) { return x * table[i].weights; });
                    if(std::abs(*std::max_element(tw2.begin(), tw2.end())) < ceps1)
                        iz2++;
                    else
                        iz2 = 0;
                } else {
                    std::fill(t2.begin(), t2.end(), 0.0);
                    std::fill(tw2.begin(), tw2.end(), 0.0);
                }

                for(size_t j = 0; j < size; ++j) {
                    sum[j] += tw1[j] + tw2[j];
                    tw1[j] = std::abs(tw1[j]);
                    tw2[j] = std::abs(tw2[j]);
                    t1[j] = std::abs(t1[j]);
                    t2[j] = std::abs(t2[j]);
                    fm1[j] = std::max({fm1[j], tw1[j], tw2[j]});
                    fm2[j] = std::max({fm2[j], t1[j], t2[j]});
                }
            }
        }

        for(size_t j = 0; j < size; ++j) s1[j] = ax * h * sum[j];
        double eps1 = *std::max_element(fm1.begin(), fm1.end()) * ceps1;
        double eps2 = *std::max_element(fm2.begin(), fm2.end()) * ceps2;
        auto diff1 = s1;
        auto diff2 = s1;
        std::transform(diff1.begin(), diff1.end(), s2.begin(), diff1.begin(), std::minus<>());
        std::transform(diff1.begin(), diff1.end(), s3.begin(), diff1.begin(), std::minus<>());
        double d1 = log10(std::abs(*std::max_element(diff1.begin(), diff1.end())));
        double d2 = log10(std::abs(*std::max_element(diff2.begin(), diff2.end())));
        double d3 = log10(eps1) - 1;
        double d4 = log10(eps2) - 1;

        if(k <= 2)
            err = 1.0;
        else if(d1 == -16)
            err = 0.0;
        else {
            double val1 = d1 * d1 / d2;
            double val2 = 2 * d1;
            double max = std::max({val1, val2, d3, d4});
            err = pow(10, static_cast<int>(std::min(0.0, max)));
        }

        spdlog::debug("Number of evaluations = {}", evals);

        if(k > 3 && err < eps1) return s1;
        if(k >= 3 && err < eps2) return s1;
    }

    return s1;
}
