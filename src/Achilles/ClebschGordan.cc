#include "Achilles/ClebschGordan.hh"
#include "Achilles/Utilities.hh"
#include "fmt/format.h"

#include <algorithm>
#include <cmath>
#include <utility>

using achilles::SpinState;

double achilles::ClebschGordan(SpinState j1, SpinState j2, SpinState j3) {
    // Ensure that momentum is conserved in z-axis
    j3.zaxis = j1.zaxis + j2.zaxis;

    if(j1.total < 0 || j2.total < 0 || j3.total < 0 || j3.total > j1.total + j2.total ||
       j3.total < std::abs(j1.total - j2.total)) {
        throw std::runtime_error(
            fmt::format("ClebschGordan: Got invalid total spins (j1={}, j2={}, j3={})", j1.total,
                        j2.total, j3.total));
    } else if(std::abs(j1.zaxis) > j1.total) {
        throw std::runtime_error(fmt::format(
            "ClebschGordan: Got invalid spin state (j = {}, m = {})", j1.total, j1.zaxis));
    } else if(std::abs(j2.zaxis) > j2.total) {
        throw std::runtime_error(fmt::format(
            "ClebschGordan: Got invalid spin state (j = {}, m = {})", j2.total, j2.zaxis));
    } else if(std::abs(j3.zaxis) > j3.total) {
        throw std::runtime_error(fmt::format(
            "ClebschGordan: Got invalid spin state (j = {}, m = {})", j3.total, j3.zaxis));
    } else if(j1.total == 0 || j2.total == 0) {
        return 1;
    }

    double sign = 1;
    if(j1.total < j2.total) {
        std::swap(j1, j2);
        sign = (j3.total - j1.total - j2.total) / 2 % 2 == 0 ? 1 : -1;
    }
    if(j3.zaxis < 0) {
        j1.zaxis = -j1.zaxis;
        j2.zaxis = -j2.zaxis;
        j3.zaxis = -j3.zaxis;
        sign = (j3.total - j1.total - j2.total) / 2 % 2 == 0 ? 1 : -1;
    }

    double numerator = (j3.total + 1) * factorial((j1.total + j2.total - j3.total) / 2) *
                       factorial((-j1.total + j2.total + j3.total) / 2) *
                       factorial((j1.total - j2.total + j3.total) / 2);
    double denominator = factorial((j1.total + j2.total + j3.total) / 2 + 1);
    double factor1 = factorial((j1.total + j1.zaxis) / 2) * factorial((j1.total - j1.zaxis) / 2);
    double factor2 = factorial((j2.total + j2.zaxis) / 2) * factorial((j2.total - j2.zaxis) / 2);
    double factor3 = factorial((j3.total + j3.zaxis) / 2) * factorial((j3.total - j3.zaxis) / 2);
    double sum_term = 0;
    int kmin =
        std::max({0, (j2.total - j3.total - j1.zaxis) / 2, (j1.total + j2.zaxis - j3.total) / 2});
    int kmax = std::min({(j2.total + j2.zaxis) / 2, (j1.total - j1.zaxis) / 2,
                         (j1.total + j2.total - j3.total) / 2});
    for(int k = kmin; k < kmax + 1; ++k) {
        auto term =
            pow(-1, k) /
            (factorial(k) * factorial((j1.total + j2.total - j3.total) / 2 - k) *
             factorial((j1.total - j1.zaxis) / 2 - k) * factorial((j2.total + j2.zaxis) / 2 - k) *
             factorial((j3.total - j2.total + j1.zaxis) / 2 + k) *
             factorial((j3.total - j1.total - j2.zaxis) / 2 + k));
        sum_term += term;
    }
    return sqrt(numerator / denominator) * sqrt(factor1 * factor2 * factor3) * sum_term * sign;
}
