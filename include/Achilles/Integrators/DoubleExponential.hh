#ifndef DOUBLEEXPONENTIAL_HH
#define DOUBLEEXPONENTIAL_HH

#include <algorithm>
#include <array>
#include <cmath>

#include "nuchic/Integrators/AdaptiveIntegrator.hh"
#include "nuchic/Utilities.hh"

namespace nuchic {
namespace Integrator {

// by Xeo, from https://stackoverflow.com/a/13294458/420683
template <size_t... Is> struct seq {};
template <size_t N, size_t... Is> struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};
template <size_t... Is> struct gen_seq<0, Is...> : seq<Is...> {};

// Generator functions
// by dyp, from https://stackoverflow.com/a/19016627/9201027
template <class Generator, size_t... Is>
constexpr auto generate_array_helper(Generator g, seq<Is...>)
    -> std::array<decltype(g(size_t{}, sizeof...(Is))), sizeof...(Is)> {
    return {{g(Is, sizeof...(Is))...}};
}

template <size_t tcount, class Generator>
constexpr auto generate_array(Generator g)
    -> decltype(generate_array_helper(g, gen_seq<tcount>{})) {
    return generate_array_helper(g, gen_seq<tcount>{});
}

// constexpr for exp using a Taylor series expansion
template <typename T> constexpr T cexp(const T &x) {
    return 1 + x + ipow<T>(x, 2) / factorial<T>(2) + ipow<T>(x, 3) / factorial<T>(3) +
           ipow<T>(x, 4) / factorial<T>(4) + ipow<T>(x, 5) / factorial<T>(5) +
           ipow<T>(x, 6) / factorial<T>(6) + ipow<T>(x, 7) / factorial<T>(7) +
           ipow<T>(x, 8) / factorial<T>(8) + ipow<T>(x, 9) / factorial<T>(9) +
           ipow<T>(x, 10) / factorial<T>(10) + ipow<T>(x, 11) / factorial<T>(11) +
           ipow<T>(x, 12) / factorial<T>(12) + ipow<T>(x, 13) / factorial<T>(13) +
           ipow<T>(x, 14) / factorial<T>(14) + ipow<T>(x, 15) / factorial<T>(15) +
           ipow<T>(x, 16) / factorial<T>(16) + ipow<T>(x, 17) / factorial<T>(17) +
           ipow<T>(x, 18) / factorial<T>(18) + ipow<T>(x, 19) / factorial<T>(19) +
           ipow<T>(x, 20) / factorial<T>(20) + ipow<T>(x, 21) / factorial<T>(21) +
           ipow<T>(x, 22) / factorial<T>(22) + ipow<T>(x, 23) / factorial<T>(23) +
           ipow<T>(x, 24) / factorial<T>(24) + ipow<T>(x, 25) / factorial<T>(25) +
           ipow<T>(x, 26) / factorial<T>(26) + ipow<T>(x, 27) / factorial<T>(27) +
           ipow<T>(x, 28) / factorial<T>(28) + ipow<T>(x, 29) / factorial<T>(29) +
           ipow<T>(x, 30) / factorial<T>(30) + ipow<T>(x, 31) / factorial<T>(31) +
           ipow<T>(x, 32) / factorial<T>(32) + ipow<T>(x, 33) / factorial<T>(33) +
           ipow<T>(x, 34) / factorial<T>(34) + ipow<T>(x, 35) / factorial<T>(35) +
           ipow<T>(x, 36) / factorial<T>(36) + ipow<T>(x, 37) / factorial<T>(37) +
           ipow<T>(x, 38) / factorial<T>(38) + ipow<T>(x, 39) / factorial<T>(39) +
           ipow<T>(x, 40) / factorial<T>(40) + ipow<T>(x, 41) / factorial<T>(41) +
           ipow<T>(x, 42) / factorial<T>(42) + ipow<T>(x, 43) / factorial<T>(43) +
           ipow<T>(x, 44) / factorial<T>(44) + ipow<T>(x, 45) / factorial<T>(45);
}

// Setup abscissa and weights for Double Exponential integration
constexpr size_t _phases = 6;
constexpr size_t _size = 6 * 1 << _phases;

struct DEPoint {
    double abscissa;
    double weights;
};

constexpr double t2(const size_t &k) {
    return cexp(static_cast<double>(k) / ipow(2, _phases));
}

constexpr double u1(const size_t &k) {
    return M_PI_4 * (t2(k) + 1.0 / t2(k));
}

constexpr double t3(const size_t &k) {
    return cexp(M_PI_4 * (t2(k) - 1.0 / t2(k)));
}

constexpr double t4(const size_t &k) {
    return (t3(k) + 1.0 / t3(k)) / 2;
}

constexpr double GetAbcs(const size_t &k) {
    return 1.0 / (t3(k) * t4(k));
}

constexpr double GetWeight(const size_t &k) {
    return u1(k) / (t4(k) * t4(k));
}

constexpr DEPoint GeneratePoint(size_t curr, size_t) {
    return {GetAbcs(curr), GetWeight(curr)};
}

// Table for Double Exponential integration
static constexpr auto table = generate_array<_size>(GeneratePoint);

class DoubleExponential : public AdaptiveIntegrator {
  public:
    DoubleExponential(bool cache = false) : AdaptiveIntegrator(cache) {}
    DoubleExponential(const FunctionS &func, bool cache = false)
        : AdaptiveIntegrator(func, cache) {}
    DoubleExponential(const FunctionV &func, bool cache = false)
        : AdaptiveIntegrator(func, cache) {}

    double Integrate(const double &, const double &, const double &, const double &) override;
    std::vector<double> IntegrateVec(const double &, const double &, const double &,
                                     const double &) override;
};

} // namespace Integrator
} // namespace nuchic

#endif
