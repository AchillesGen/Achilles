#include "nuchic/Interpolation.hh"

#include <iostream>

extern "C" {
    // Fortran 1D interpolation interface
    nuchic::Interp1D *CreateInterp1D(const double *x, const double *y, int n) {
        auto *interp = new nuchic::Interp1D(); 

        std::vector<double> xVec(x, x+n);
        std::vector<double> yVec(y, y+n);

        interp -> CubicSpline(xVec, yVec);
        return interp;
    }

    void DeleteInterp1d(nuchic::Interp1D *interp) {
        delete interp;
    }

    double Interp1DMin(nuchic::Interp1D *interp) {
        return interp -> min();
    }

    double Interp1DMax(nuchic::Interp1D *interp) {
        return interp -> max();
    }

    double Interpolate1D(nuchic::Interp1D *interp, const double x) {
        return interp -> operator()(x); 
    }

    // Fortran 2D interpolation interface
    nuchic::Interp2D *CreateInterp2D(double *x, double *y, double *z, int n1, int n2) {
        auto *interp = new nuchic::Interp2D(); 

        std::vector<double> xVec(x, x+n1);
        std::vector<double> yVec(y, y+n2);
        std::vector<double> zVec(z, z+n1*n2);

        interp -> BicubicSpline(xVec, yVec, zVec);
        return interp;
    }

    void DeleteInterp2D(nuchic::Interp2D *interp) {
        delete interp;
    }

    double Interp2DXMin(nuchic::Interp2D *interp) {
        return interp -> Xmin();
    }

    double Interp2DXMax(nuchic::Interp2D *interp) {
        return interp -> Xmax();
    }

    double Interp2DYMin(nuchic::Interp2D *interp) {
        return interp -> Ymin();
    }

    double Interp2DYMax(nuchic::Interp2D *interp) {
        return interp -> Ymax();
    }

    double Interpolate2D(nuchic::Interp2D *interp, double x, double y) {
        return interp -> operator()(x, y);
    }
}
