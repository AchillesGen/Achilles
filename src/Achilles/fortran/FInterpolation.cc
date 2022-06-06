#include "fmt/ostream.h"
#include "Achilles/Interpolation.hh"

#include <iostream>

extern "C" {
    // Fortran 1D interpolation interface
    achilles::Interp1D *CreateInterp1D(const double *x, const double *y, int n, int mode) {
        std::vector<double> xVec(x, x+n);
        std::vector<double> yVec(y, y+n);

        achilles::Interp1D *interp = nullptr;
        switch(mode) {
            case 0:
                interp = new achilles::Interp1D(xVec, yVec,
                        achilles::InterpolationType::NearestNeighbor); 
                break;
            case 1:
            case 2:
            case 3:
                interp = new achilles::Interp1D(xVec, yVec,
                        achilles::InterpolationType::Polynomial); 
                interp -> SetPolyOrder(static_cast<size_t>(mode));
                break;
            case 4:
                interp = new achilles::Interp1D(xVec, yVec,
                        achilles::InterpolationType::CubicSpline); 
                interp -> CubicSpline();
                break;
            default:
                throw std::runtime_error("Invalid interpolation mode");
        }

        return interp;
    }

    void DeleteInterp1d(achilles::Interp1D *interp) {
        delete interp;
    }

    double Interp1DMin(achilles::Interp1D *interp) {
        return interp -> min();
    }

    double Interp1DMax(achilles::Interp1D *interp) {
        return interp -> max();
    }

    double Interpolate1D(achilles::Interp1D *interp, const double x) {
        return interp -> operator()(x); 
    }

    // Fortran 2D interpolation interface
    achilles::Interp2D *CreateInterp2D(double *x, double *y, double *z, int n1, int n2, int mode) {
        std::vector<double> xVec(x, x+n1);
        std::vector<double> yVec(y, y+n2);
        std::vector<double> zVec(z, z+n1*n2);

        achilles::Interp2D *interp;
        switch(mode) {
            case 0:
                interp = new achilles::Interp2D(xVec, yVec, zVec,
                        achilles::InterpolationType::NearestNeighbor); 
                break;
            case 1:
            case 2:
            case 3:
                interp = new achilles::Interp2D(xVec, yVec, zVec,
                        achilles::InterpolationType::Polynomial); 
                interp -> SetPolyOrder(static_cast<size_t>(mode), static_cast<size_t>(mode));
                break;
            case 4:
                interp = new achilles::Interp2D(xVec, yVec, zVec,
                        achilles::InterpolationType::CubicSpline); 
                interp -> BicubicSpline();
                break;
            default:
                throw std::runtime_error("Invalid interpolation mode");
        }

        return interp;
    }

    void DeleteInterp2D(achilles::Interp2D *interp) {
        delete interp;
    }

    double Interp2DXMin(achilles::Interp2D *interp) {
        return interp -> Xmin();
    }

    double Interp2DXMax(achilles::Interp2D *interp) {
        return interp -> Xmax();
    }

    double Interp2DYMin(achilles::Interp2D *interp) {
        return interp -> Ymin();
    }

    double Interp2DYMax(achilles::Interp2D *interp) {
        return interp -> Ymax();
    }

    double Interpolate2D(achilles::Interp2D *interp, double x, double y) {
        try {
            return interp -> operator()(x, y);
        } catch (std::domain_error &e) {
            return 0;
        }
    }
}
