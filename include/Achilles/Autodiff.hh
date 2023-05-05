#ifndef AUTODIFF_HH
#define AUTODIFF_HH

#include <cmath>
#include <type_traits>
#include <utility>

namespace achilles {

using std::abs;
using std::cos;
using std::cosh;
using std::exp;
using std::log;
using std::pow;
using std::sin;
using std::sinh;
using std::tan;

class Dual {
  private:
    double m_val{};
    double m_eps{};

  public:
    Dual() = default;
    explicit Dual(double val) : m_val{val}, m_eps{1} {}
    Dual(double val, double eps) : m_val{val}, m_eps{eps} {}

    // Getters and Setters
    double Value() const { return m_val; }
    double &Value() { return m_val; }
    double Derivative() const { return m_eps; }
    double &Derivative() { return m_eps; }
    std::pair<double, double> Get() { return {m_val, m_eps}; }
    explicit operator double() const { return m_val; }
};

// Operators
// Addition:
Dual operator+(const Dual &, const Dual &);
Dual operator+(const Dual &);

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Dual operator+(const Dual &x, const T &y) {
    Dual result = x;
    result.Value() += y;
    return result;
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Dual operator+(const T &x, const Dual &y) {
    Dual result = y;
    result.Value() += x;
    return result;
}

// Subtraction:
Dual operator-(const Dual &, const Dual &);
Dual operator-(const Dual &);

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Dual operator-(const Dual &x, const T &y) {
    return {x.Value() - y, x.Derivative()};
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Dual operator-(const T &x, const Dual &y) {
    return {x - y.Value(), -y.Derivative()};
}

// Multiplication:
Dual operator*(const Dual &, const Dual &);

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Dual operator*(const Dual &x, const T &y) {
    return {x.Value() * y, x.Derivative() * y};
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Dual operator*(const T &x, const Dual &y) {
    return {y.Value() * x, y.Derivative() * x};
}

// Division:
Dual operator/(const Dual &, const Dual &);

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Dual operator/(const Dual &x, const T &y) {
    return {x.Value() / y, x.Derivative() / y};
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Dual operator/(const T &x, const Dual &y) {
    return {x / y.Value(), -x * y.Derivative() / std::pow(y.Value(), 2)};
}

// Other functions
Dual sin(const Dual &);
Dual cos(const Dual &);
Dual tan(const Dual &);
Dual exp(const Dual &);
Dual log(const Dual &);
Dual abs(const Dual &);
Dual cosh(const Dual &);
Dual sinh(const Dual &);
Dual sech(const Dual &);

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Dual pow(const Dual &x, const T &power) {
    return {std::pow(x.Value(), power), power * std::pow(x.Value(), power - 1) * x.Derivative()};
}

} // namespace achilles

#endif
