#include "fmt/ostream.h"
#include "nuchic/Interpolation.hh"

extern "C" {
    // Fortran 1D interpolation interface
    nuchic::Interp1D *CreateInterp1D(const double *x, const double *y, int n, int mode) {
        std::vector<double> xVec(x, x+n);
        std::vector<double> yVec(y, y+n);

        nuchic::Interp1D *interp = nullptr;
        switch(mode) {
            case 0:
                interp = new nuchic::Interp1D(xVec, yVec,
                        nuchic::InterpolationType::NearestNeighbor); 
                break;
            case 1:
            case 2:
            case 3:
                interp = new nuchic::Interp1D(xVec, yVec,
                        nuchic::InterpolationType::Polynomial); 
                interp -> SetPolyOrder(static_cast<size_t>(mode));
                break;
            case 4:
                interp = new nuchic::Interp1D(xVec, yVec,
                        nuchic::InterpolationType::CubicSpline); 
                interp -> CubicSpline();
                break;
            default:
                throw std::runtime_error("Invalid interpolation mode");
        }

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
    nuchic::Interp2D *CreateInterp2D(double *x, double *y, double *z, int n1, int n2, int mode) {
        std::vector<double> xVec(x, x+n1);
        std::vector<double> yVec(y, y+n2);
        std::vector<double> zVec(z, z+n1*n2);

        nuchic::Interp2D *interp;
        switch(mode) {
            case 0:
                interp = new nuchic::Interp2D(xVec, yVec, zVec,
                        nuchic::InterpolationType::NearestNeighbor); 
                break;
            case 1:
            case 2:
            case 3:
                interp = new nuchic::Interp2D(xVec, yVec, zVec,
                        nuchic::InterpolationType::Polynomial); 
                interp -> SetPolyOrder(static_cast<size_t>(mode), static_cast<size_t>(mode));
                break;
            case 4:
                interp = new nuchic::Interp2D(xVec, yVec, zVec,
                        nuchic::InterpolationType::CubicSpline); 
                interp -> BicubicSpline();
                break;
            default:
                throw std::runtime_error("Invalid interpolation mode");
        }

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
        try {
            return interp -> operator()(x, y);
        } catch (std::domain_error &e) {
            return 0;
        }
    }
}
