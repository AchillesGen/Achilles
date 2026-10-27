// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <iostream>

#include "Achilles/Interpolation.hh"
#include "fmt/format.h"

constexpr double ipow(double x, size_t exponent) {
    return (exponent == 0)       ? 1
           : (exponent % 2 == 0) ? ipow(x * x, exponent / 2)
                                 : x * ipow(x * x, (exponent - 1) / 2);
}

using namespace achilles;

double achilles::NevilleInterpolate(const std::vector<double> &nodes,
                                    const std::vector<double> &values, size_t count, double x) {
    // Neville's algorithm (E. H. Neville, J. Indian Math. Soc. 20 (1934) 87): evaluate the unique
    // degree-(count-1) polynomial through the given nodes at x, by repeatedly combining adjacent
    // lower-degree interpolants using
    //
    //     P_{i,i}   = values[i]
    //     P_{i,i+k} = [ (x - nodes[i+k]) P_{i,i+k-1} + (nodes[i] - x) P_{i+1,i+k} ]
    //                 / (nodes[i] - nodes[i+k])
    //
    // The tableau is carried in one rolling buffer: after pass k, tableau[i] holds P_{i,i+k}.
    if(count == 0) throw std::runtime_error("NevilleInterpolate: no interpolation nodes");
    if(values.size() < count || nodes.size() < count)
        throw std::runtime_error("NevilleInterpolate: fewer nodes than requested");

    std::vector<double> tableau(values.begin(),
                                values.begin() + static_cast<std::ptrdiff_t>(count));
    for(size_t k = 1; k < count; ++k) {
        for(size_t i = 0; i + k < count; ++i) {
            const double spread = nodes[i] - nodes[i + k];
            if(spread == 0)
                throw std::runtime_error("NevilleInterpolate: repeated interpolation node");
            tableau[i] =
                ((x - nodes[i + k]) * tableau[i] + (nodes[i] - x) * tableau[i + 1]) / spread;
        }
    }
    return tableau[0];
}

constexpr double Interp1D::maxDeriv;

Interp1D::Interp1D(const std::vector<double> &x, const std::vector<double> &y,
                   InterpolationType mode)
    : kMode{mode} {
    if(!std::is_sorted(x.begin(), x.end())) throw std::runtime_error("Inputs must be increasing.");
    if(std::adjacent_find(x.begin(), x.end()) != x.end())
        throw std::runtime_error("Inputs must all be unique.");
    if(x.size() != y.size())
        throw std::runtime_error("Input and output arrays must be the same size.");

    knotX = x;
    knotY = y;
}

void Interp1D::CubicSpline(const double &derivLeft, const double &derivRight) {
    // The natural/clamped interpolating cubic spline is fixed by its second derivatives M_i,
    // which satisfy the tridiagonal system
    //
    //     h_{i-1} M_{i-1} + 2 (h_{i-1} + h_i) M_i + h_i M_{i+1}
    //         = 6 [ (y_{i+1} - y_i) / h_i - (y_i - y_{i-1}) / h_{i-1} ],    i = 1 .. n-2
    //
    // with h_i = x_{i+1} - x_i, closed at each end by either a natural condition (M = 0) or a
    // clamped one (prescribed first derivative).  Solved here by the Thomas algorithm.
    const std::size_t n = knotX.size();
    derivs2.assign(n, 0.0);
    if(n < 2) {
        kSplineInit = true;
        return;
    }

    const auto h = [&](std::size_t i) { return knotX[i + 1] - knotX[i]; };
    const auto slope = [&](std::size_t i) { return (knotY[i + 1] - knotY[i]) / h(i); };

    std::vector<double> sub(n, 0.0), diag(n, 0.0), sup(n, 0.0), rhs(n, 0.0);

    if(derivLeft >= maxDeriv) { // natural
        diag[0] = 1.0;
    } else { // clamped
        diag[0] = 2.0 * h(0);
        sup[0] = h(0);
        rhs[0] = 6.0 * (slope(0) - derivLeft);
    }

    for(std::size_t i = 1; i + 1 < n; ++i) {
        sub[i] = h(i - 1);
        diag[i] = 2.0 * (h(i - 1) + h(i));
        sup[i] = h(i);
        rhs[i] = 6.0 * (slope(i) - slope(i - 1));
    }

    if(derivRight >= maxDeriv) { // natural
        diag[n - 1] = 1.0;
    } else { // clamped
        sub[n - 1] = h(n - 2);
        diag[n - 1] = 2.0 * h(n - 2);
        rhs[n - 1] = 6.0 * (derivRight - slope(n - 2));
    }

    for(std::size_t i = 1; i < n; ++i) { // forward elimination
        const double factor = sub[i] / diag[i - 1];
        diag[i] -= factor * sup[i - 1];
        rhs[i] -= factor * rhs[i - 1];
    }
    derivs2[n - 1] = rhs[n - 1] / diag[n - 1]; // back substitution
    for(std::size_t i = n - 1; i > 0; --i)
        derivs2[i - 1] = (rhs[i - 1] - sup[i - 1] * derivs2[i]) / diag[i - 1];

    kSplineInit = true;
}

double Interp1D::operator()(const double &x) const {
    // Ensure the interpolation is initialized first
    if(!kSplineInit && kMode == InterpolationType::CubicSpline)
        throw std::runtime_error("Interpolation is not initialized!");

    // Disallow extrapolation
    if(x > knotX.back())
        throw std::domain_error(
            fmt::format("Input ({}) greater than maximum value ({})", x, knotX.back()));
    if(x < knotX.front())
        throw std::domain_error(
            fmt::format("Input ({}) less than minimum value ({})", x, knotX.front()));

    // Find range by binary_search
    auto idxHigh = static_cast<size_t>(
        std::distance(knotX.begin(), std::upper_bound(knotX.begin(), knotX.end(), x)));
    auto idxLow = idxHigh - 1;

    double result = 0;
    switch(kMode) {
    case InterpolationType::NearestNeighbor:
        result = x - knotX[idxLow] < knotX[idxHigh] - x ? knotY[idxLow] : knotY[idxHigh];
        break;
    case InterpolationType::Polynomial:
        result = PolynomialInterp(x);
        break;
    case InterpolationType::CubicSpline:
        const double height = knotX[idxHigh] - knotX[idxLow];
        const double a = (knotX[idxHigh] - x) / height;
        const double b = (x - knotX[idxLow]) / height;

        result = a * knotY[idxLow] + b * knotY[idxHigh] +
                 ((ipow(a, 3) - a) * derivs2[idxLow] + (ipow(b, 3) - b) * derivs2[idxHigh]) *
                     ipow(height, 2) / 6.0;
        break;
    }
    return result;
}

// double Interp1D::Integrate(double a, double b) const {
//     if(a < knotX.front() || a > knotX.back() || b < knotX.front() || b > knotX.back())
//         throw std::domain_error(
//             fmt::format("Invalid integration region [{}, {}] for function defined on [{}, {}]",
//                         a, b, knotX.front(), knotX.back()));
//     if(b < a) return -Integrate(b, a);
// }

double Interp1D::Integrate() const {
    double result = 0;
    switch(kMode) {
    case InterpolationType::NearestNeighbor:
        for(size_t i = 0; i < knotX.size() - 1; ++i) {
            double width = knotX[i + 1] - knotX[i];
            result += width / 2 * (knotY[i] + knotY[i + 1]);
        }
        break;
    case InterpolationType::Polynomial:
        if(polyOrder == 2) {
            for(size_t i = 0; i < knotX.size() - 1; ++i) {
                double width = knotX[i + 1] - knotX[i];
                result += width / 2 * (knotY[i] + knotY[i + 1]);
            }
        } else if(polyOrder == 3) {
            for(size_t i = 0; i < knotX.size() - 2; i += 2) {
                double h1 = knotX[i + 1] - knotX[i];
                double h2 = knotX[i + 2] - knotX[i + 1];
                double f0 = knotY[i];
                double f1 = knotY[i + 1];
                double f2 = knotY[i + 2];
                result +=
                    (h1 + h2) *
                    (-f2 * h1 * (h1 - 2 * h2) + f0 * (2 * h1 - h2) * h2 + f1 * pow(h1 + h2, 2)) /
                    (6 * h1 * h2);
            }
        } else {
            throw std::runtime_error("Not Implemented!");
        }
        break;
    case InterpolationType::CubicSpline:
        for(size_t i = 0; i < knotX.size() - 1; ++i) {
            result += -(knotX[i + 1] - knotX[i]) *
                      ((derivs2[i + 1] + derivs2[i]) * pow(knotX[i + 1] - knotX[i], 2) -
                       12 * (knotY[i + 1] + knotY[i])) /
                      24.0;
        }
        break;
    }
    return result;
}

double Interp1D::PolynomialInterp(double x) const {
    auto idx = static_cast<size_t>(
        std::distance(knotX.begin(), std::lower_bound(knotX.begin(), knotX.end(), x)));
    while(idx < polyOrder / 2) ++idx;
    while(knotX.size() - idx < polyOrder / 2 + polyOrder % 2) --idx;
    std::vector<double> xInterp(&knotX[idx] - polyOrder / 2,
                                &knotX[idx] + polyOrder / 2 + polyOrder % 2);
    std::vector<double> tmp(&knotY[idx] - polyOrder / 2,
                            &knotY[idx] + polyOrder / 2 + polyOrder % 2);

    return NevilleInterpolate(xInterp, tmp, static_cast<size_t>(polyOrder), x);
}

Interp2D::Interp2D(const std::vector<double> &x, const std::vector<double> &y,
                   const std::vector<double> &z, InterpolationType mode)
    : kMode{mode} {
    if(!std::is_sorted(x.begin(), x.end())) throw std::runtime_error("Inputs must be increasing.");
    if(std::adjacent_find(x.begin(), x.end()) != x.end())
        throw std::runtime_error("Inputs must all be unique.");
    if(!std::is_sorted(y.begin(), y.end())) throw std::runtime_error("Inputs must be increasing.");
    if(std::adjacent_find(y.begin(), y.end()) != y.end())
        throw std::runtime_error("Inputs must all be unique.");
    if(x.size() * y.size() != z.size())
        throw std::runtime_error("Input and output arrays must be the same size.");

    knotX = x;
    knotY = y;
    knotZ = z;
}

void Interp2D::BicubicSpline() {
    for(std::size_t i = 0; i < knotX.size(); ++i) {
        derivs2.emplace_back(
            knotY, std::vector<double>(knotZ.begin() + static_cast<int>(i * knotY.size()),
                                       knotZ.begin() + static_cast<int>((i + 1) * knotY.size())));
        derivs2.back().CubicSpline();
    }

    kSplineInit = true;
}

double Interp2D::operator()(const double &x, const double &y) const {
    // Ensure the interpolation is initialized first
    if(!kSplineInit && kMode == InterpolationType::CubicSpline)
        throw std::runtime_error("Interpolation is not initialized!");

    // Disallow extrapolation
    if(x > knotX.back())
        throw std::domain_error(
            fmt::format("Input ({}) greater than maximum x value ({})", x, knotX.back()));
    if(x < knotX.front())
        throw std::domain_error(
            fmt::format("Input ({}) less than minimum x value ({})", x, knotX.front()));
    if(y > knotY.back())
        throw std::domain_error(
            fmt::format("Input ({}) greater than maximum y value ({})", y, knotY.back()));
    if(y < knotY.front())
        throw std::domain_error(
            fmt::format("Input ({}) less than minimum y value ({})", y, knotY.front()));

    double result = 0;
    switch(kMode) {
    case InterpolationType::NearestNeighbor:
        result = NearestNeighbor(x, y);
        break;
    case InterpolationType::Polynomial:
        result = PolynomialInterp(x, y);
        break;
    case InterpolationType::CubicSpline:
        std::vector<double> zTmp(knotX.size());
        for(std::size_t i = 0; i < knotX.size(); ++i) zTmp[i] = derivs2[i](y);

        Interp1D interp(knotX, zTmp);
        interp.SetType(InterpolationType::CubicSpline);
        interp.CubicSpline();
        result = interp(x);
        break;
    }

    return result;
}

double Interp2D::NearestNeighbor(double x, double y) const {
    // Find range by binary_search
    auto idxHighX = static_cast<size_t>(
        std::distance(knotX.begin(), std::upper_bound(knotX.begin(), knotX.end(), x)));
    auto idxLowX = idxHighX - 1;
    auto idxHighY = static_cast<size_t>(
        std::distance(knotY.begin(), std::upper_bound(knotY.begin(), knotY.end(), y)));
    auto idxLowY = idxHighY - 1;
    auto idxX = x - knotX[idxLowX] < knotX[idxHighX] - x ? idxLowX : idxHighX;
    auto idxY = y - knotY[idxLowY] < knotY[idxHighY] - y ? idxLowY : idxHighY;
    return knotZ[idxY + knotY.size() * idxX];
}

double Interp2D::PolynomialInterp(double x, double y) const {
    // Find point in x direction
    auto idxX = static_cast<size_t>(
        std::distance(knotX.begin(), std::lower_bound(knotX.begin(), knotX.end(), x)));
    while(idxX < polyOrderX / 2) ++idxX;
    while(knotX.size() - idxX < polyOrderX / 2 + polyOrderX % 2) --idxX;
    std::vector<double> xInterp(&knotX[idxX] - polyOrderX / 2,
                                &knotX[idxX] + polyOrderX / 2 + polyOrderX % 2);

    // Find point in y direction
    auto idxY = static_cast<size_t>(
        std::distance(knotY.begin(), std::lower_bound(knotY.begin(), knotY.end(), y)));
    while(idxY < polyOrderY / 2) ++idxY;
    while(knotY.size() - idxY < polyOrderY / 2 + polyOrderY % 2) --idxY;
    std::vector<double> yInterp(&knotY[idxY] - polyOrderY / 2,
                                &knotY[idxY] + polyOrderY / 2 + polyOrderY % 2);

    std::vector<double> tmp(polyOrderY);
    std::vector<double> tmp2(polyOrderX);
    for(size_t i = 0; i < polyOrderX; ++i) {
        for(size_t j = 0; j < polyOrderY; ++j) {
            tmp[j] = knotZ[idxY + j - polyOrderY / 2 + knotY.size() * (idxX - polyOrderX / 2 + i)];
        }
        tmp2[i] = NevilleInterpolate(yInterp, tmp, polyOrderY, y);
    }
    return NevilleInterpolate(xInterp, tmp2, polyOrderX, x);
}
