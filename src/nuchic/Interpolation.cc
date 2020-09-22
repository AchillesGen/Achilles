#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <iostream>

#include "fmt/format.h"
#include "nuchic/Interpolation.hh"

using namespace nuchic;

constexpr double Interp1D::maxDeriv;
void Interp1D::CubicSpline(const std::vector<double>& x, const std::vector<double>& y,
                                   const double& derivLeft, const double& derivRight) {
    
    if(!std::is_sorted(x.begin(), x.end()))
        std::runtime_error("Inputs must be increasing.");
    if(std::adjacent_find(x.begin(), x.end()) != x.end())
        std::runtime_error("Inputs must all be unique.");
    if(x.size() != y.size())
        std::runtime_error("Input and output arrays must be the same size.");

    knotX = x;
    knotY = y;
    const std::size_t n = knotX.size();
    std::vector<double> u(n);
    derivs2.resize(n);

    if(derivLeft >= maxDeriv) {
        derivs2[0] = 0.0;
        u[0] = 0.0;
    } else {
        derivs2[0] = -0.5;
        u[0] = (3/(knotX[1]-knotX[0]))*((knotY[1]-knotY[0])/(knotX[1]-knotX[0])-derivLeft);
    }

    for(std::size_t i = 1; i < n-1; ++i) {
        double sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
        double p = sig*derivs2[i-1]+2;
        derivs2[i] = (sig-1.0)/p;
        u[i] = (knotY[i+1]-knotY[i])/(knotX[i+1]-knotX[i])
            - (knotY[i]-knotY[i-1])/(knotX[i]-knotX[i-1]);
        u[i] = (6.0*u[i]/(knotX[i+1]-knotX[i-1])-sig*u[i-1])/p;
    }

    double dn{}, un{};
    if(derivRight >= maxDeriv) {
        dn = 0.0;
        un = 0.0;
    } else {
        dn = 0.5;
        un = (3/(knotX[n-1]-knotX[n-2]))*(derivRight-(knotY[n-1]-knotY[n-2])
                /(knotX[n-1]-knotX[n-2]));
    }

    derivs2[n-1]=(un-dn*u[n-2])/(dn*derivs2[n-2]+1.0);
    for(std::size_t i = n-1; i > 0; --i) {
        derivs2[i-1] = derivs2[i-1]*derivs2[i]+u[i-1];
    }

    kInit = true;
}

double Interp1D::operator()(const double& x) const {
    // Ensure the interpolation is initialized first
    if(!kInit)
        throw std::runtime_error("Interpolation is not initialized!");

    // Disallow extrapolation
    if(x > knotX.back()) 
        throw std::domain_error(fmt::format("Input ({}) greater than maximum value ({})", x, knotX.back()));

    std::size_t idxLow = 0, idxHigh = knotX.size(), idx = 0; 

    // Find range by bisection
    while(idxHigh - idxLow > 1) {
        idx = (idxHigh + idxLow) >> 1;
        if(knotX[idx] > x) idxHigh = idx;
        else idxLow = idx;
    }

    const double height = knotX[idxHigh] - knotX[idxLow];
    const double a = (knotX[idxHigh] - x)/height;
    const double b = (x - knotX[idxLow])/height;

    return a*knotY[idxLow] + b*knotY[idxHigh] + ((pow(a, 3) - a)*derivs2[idxLow]
            + (pow(b, 3) - b)*derivs2[idxHigh])*pow(height, 2)/6.0;
}

void Interp2D::BicubicSpline(const std::vector<double>& x, const std::vector<double>& y,
                                     const std::vector<double>& z) {
    if(!std::is_sorted(x.begin(), x.end()))
        std::runtime_error("Inputs must be increasing.");
    if(std::adjacent_find(x.begin(), x.end()) != x.end())
        std::runtime_error("Inputs must all be unique.");
    if(!std::is_sorted(y.begin(), y.end()))
        std::runtime_error("Inputs must be increasing.");
    if(std::adjacent_find(y.begin(), y.end()) != y.end())
        std::runtime_error("Inputs must all be unique.");
    if(x.size()*y.size() != z.size())
        std::runtime_error("Input and output arrays must be the same size.");

    knotX = x;
    knotY = y;
    knotZ = z;
    derivs2.resize(y.size());

    for(std::size_t i = 0; i < x.size(); ++i) {
        derivs2.emplace_back();
        derivs2[i].CubicSpline(y, std::vector<double>(z.begin()+static_cast<int>(i*y.size()), z.begin()+static_cast<int>((i+1)*y.size())));
    }

    kInit = true;
}

double Interp2D::operator()(const double& x, const double& y) const {
    // Ensure the interpolation is initialized first
    if(!kInit)
        throw std::runtime_error("Interpolation is not initialized!");

    // Disallow extrapolation
    if(x > knotX.back()) 
        throw std::domain_error(fmt::format("Input ({}) greater than maximum x value ({})", x, knotX.back()));
    if(y > knotY.back()) 
        throw std::domain_error(fmt::format("Input ({}) greater than maximum y value ({})", y, knotY.back()));


    std::vector<double> zTmp(knotX.size()); 
    for(std::size_t i = 0; i < knotX.size(); ++i)
        zTmp[i] = derivs2[i](y);

    Interp1D interp;
    interp.CubicSpline(knotX, zTmp);

    return interp(x);
}
