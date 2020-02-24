#ifndef UNITS_HH
#define UNITS_HH

namespace nuchic {
    // Time literals
    constexpr double operator"" _s(long double x) {
        return x;
    }
    constexpr double operator"" _s(unsigned long long int x) {
        return x;
    }

    // Distance literals
    constexpr double operator"" _fm(long double x) {
        return x;
    }
    constexpr double operator"" _fm(unsigned long long int x) {
        return x;
    }
    constexpr double operator"" _m(long double x) {
        return x*1.0e15;
    }
    constexpr double operator"" _m(unsigned long long int x) {
        return x*1.0e15;
    }

    // Cross-section literals
    constexpr double operator"" _mb(long double x) {
        return x;
    }
    constexpr double operator"" _mb(unsigned long long int x) {
        return x;
    }

    // Energy literals
    constexpr double operator"" _MeV(long double x) {
        return x;
    }
    constexpr double operator"" _MeV(unsigned long long int x) {
        return x;
    }
    constexpr double operator"" _GeV(long double x) {
        return x*1000;
    }
    constexpr double operator"" _GeV(unsigned long long int x) {
        return x*1000;
    }

}

#endif
