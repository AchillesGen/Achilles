#include <Achilles/Autodiff.hh>

using achilles::Dual;

Dual achilles::operator+(const Dual &x, const Dual &y) {
    return {x.Value() + y.Value(), x.Derivative() + y.Derivative()};
}

Dual achilles::operator+(const Dual &x) {
    return x;
}

Dual achilles::operator-(const Dual &x, const Dual &y) {
    return {x.Value() - y.Value(), x.Derivative() - y.Derivative()};
}

Dual achilles::operator-(const Dual &x) {
    return {-x.Value(), -x.Derivative()};
}

Dual achilles::operator*(const Dual &x, const Dual &y) {
    return {x.Value() * y.Value(), y.Value() * x.Derivative() + x.Value() * y.Derivative()};
}

Dual achilles::operator/(const Dual &x, const Dual &y) {
    return {x.Value() / y.Value(),
            (y.Value() * x.Derivative() - x.Value() * y.Derivative()) / y.Value() / y.Value()};
}

Dual achilles::sin(const Dual &x) {
    return {std::sin(x.Value()), std::cos(x.Value()) * x.Derivative()};
}

Dual achilles::cos(const Dual &x) {
    return {std::cos(x.Value()), -std::sin(x.Value()) * x.Derivative()};
}

Dual achilles::tan(const Dual &x) {
    return {std::tan(x.Value()), 1.0 / std::pow(std::cos(x.Value()), 2) * x.Derivative()};
}

Dual achilles::exp(const Dual &x) {
    return {std::exp(x.Value()), std::exp(x.Value()) * x.Derivative()};
}

Dual achilles::log(const Dual &x) {
    return {std::log(x.Value()), x.Derivative() / x.Value()};
}

Dual achilles::abs(const Dual &x) {
    double sign = x.Value() == 0 ? 0 : x.Value() / std::abs(x.Value());
    return {std::abs(x.Value()), x.Derivative() * sign};
}

Dual achilles::cosh(const Dual &x) {
    return {std::cosh(x.Value()), x.Derivative() * std::sinh(x.Value())};
}

Dual achilles::sinh(const Dual &x) {
    return {std::sinh(x.Value()), x.Derivative() * std::cosh(x.Value())};
}

Dual achilles::sech(const Dual &x) {
    return 1.0 / cosh(x);
}
