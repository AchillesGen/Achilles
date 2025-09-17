#ifndef ACHILLES_DISTRIBUTIONS
#define ACHILLES_DISTRIBUTIONS

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace achilles::distributions {

inline double BlattWeisskopf(double x, size_t angular_mom) {
    switch(angular_mom) {
    case 0:
        return 1;
    case 1:
        return x / sqrt(1 + x * x);
    case 2:
        return x * x / sqrt(9 + 3 * x * x + pow(x, 4));
    case 3:
        return x * x * x / sqrt(225 + 45 * x * x + 6 * pow(x, 4) + pow(x, 6));
    case 4:
        return pow(x, 4) /
               sqrt(11025 + 1575 * x * x + 135 * pow(x, 4) + 10 * pow(x, 6) + pow(x, 8));
    default:
        throw std::runtime_error("BlattWeisskopf: Only implemented up to angular momentum of 4");
    }
}

} // namespace achilles::distributions

#endif
