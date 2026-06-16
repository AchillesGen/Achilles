#ifndef UNITS_HH
#define UNITS_HH

#include <cmath>
#include <stdexcept>
#include <string>

#include "spdlog/spdlog.h"

namespace achilles {
constexpr double base_to_femto = 1.0e15;

// Time literals
constexpr double operator"" _s(long double x) {
    return static_cast<double>(x);
}
constexpr double operator"" _s(unsigned long long int x) {
    return static_cast<double>(x);
}

// Distance literals
constexpr double operator"" _fm(long double x) {
    return static_cast<double>(x);
}
constexpr double operator"" _fm(unsigned long long int x) {
    return static_cast<double>(x);
}
constexpr double operator"" _m(long double x) {
    return static_cast<double>(x) * base_to_femto;
}
constexpr double operator"" _m(unsigned long long int x) {
    return static_cast<double>(x) * base_to_femto;
}

// Cross-section literals
constexpr double operator"" _mb(long double x) {
    return static_cast<double>(x);
}
constexpr double operator"" _mb(unsigned long long int x) {
    return static_cast<double>(x);
}

// Energy literals
constexpr double operator"" _MeV(long double x) {
    return static_cast<double>(x);
}
constexpr double operator"" _MeV(unsigned long long int x) {
    return static_cast<double>(x);
}
constexpr double operator"" _GeV(long double x) {
    return static_cast<double>(x * 1000);
}
constexpr double operator"" _GeV(unsigned long long int x) {
    return static_cast<double>(x) * 1000;
}
constexpr double operator"" _J(long double x) {
    return static_cast<double>(x) / 1.602176634e-13;
}

// Angle literals
constexpr double operator"" _rad(long double x) {
    return static_cast<double>(x);
}
constexpr double operator"" _deg(long double x) {
    constexpr double ToRads = M_PI / 180;
    return static_cast<double>(x) * ToRads;
}

// Display units for cross-section / integral results. LAST is a sentinel marking the
// end of the valid range so the units can be iterated in one place.
enum class UnitsEnum { mub, nb, pb, em38cm2, LAST };

inline std::string ToString(UnitsEnum unit) {
    switch(unit) {
    case UnitsEnum::mub:
        return "mub";
    case UnitsEnum::nb:
        return "nb";
    case UnitsEnum::pb:
        return "pb";
    case UnitsEnum::em38cm2:
        return "e-38cm2";
    default:
        throw std::runtime_error("Invalid UnitsEnum value");
    }
}

inline double UnitScale(UnitsEnum unit) {
    switch(unit) {
    case UnitsEnum::mub:
        return 1e-3;
    case UnitsEnum::nb:
        return 1;
    case UnitsEnum::pb:
        return 1e3;
    case UnitsEnum::em38cm2:
        return 1e5;
    default:
        throw std::runtime_error("Invalid UnitsEnum value");
    }
}

// Maps a run-card string to a display unit, defaulting to nb (with a warning) when the
// string is not recognized. This is the only place that produces the default unit; every
// valid UnitsEnum is otherwise handled explicitly above (invalid values throw).
inline UnitsEnum ToUnitEnum(const std::string &unit) {
    for(int i = 0; i < static_cast<int>(UnitsEnum::LAST); ++i) {
        auto candidate = static_cast<UnitsEnum>(i);
        if(unit == ToString(candidate)) return candidate;
    }

    spdlog::warn("Unit {} not recognized. Available display units include:", unit);
    for(int i = 0; i < static_cast<int>(UnitsEnum::LAST); ++i)
        spdlog::warn("  {}", ToString(static_cast<UnitsEnum>(i)));
    spdlog::warn("Defaulting back to nb");
    return UnitsEnum::nb;
}

} // namespace achilles

#endif
