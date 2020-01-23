#ifndef INTERPOLATION_HH
#define INTERPOLATION_HH

#include <vector>

#include "pybind11/numpy.h"

namespace py = pybind11;

using pyArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

class Interp1D {
    public:
        Interp1D() : kInit(false) {}
        ~Interp1D() {}

        // Functions
        void CubicSpline(const std::vector<double>&, const std::vector<double>&,
                const double& derivLeft=1.e30, const double& derivRight=1.e30);

        double operator()(const double&) const;

    private:
        bool kInit;
        std::vector<double> knotX, knotY, derivs2;
};

class Interp2D {
    public:
        Interp2D() : kInit(false) {}
        ~Interp2D() {}

        // Functions
        void BicubicSpline(const std::vector<double>&, const std::vector<double>&,
                const std::vector<double>&);
        void BicubicSpline(const std::vector<double>&, const std::vector<double>&,
                const pyArray&);

        double operator()(const double&, const double&) const;

    private:
        bool kInit;
        std::vector<double> knotX, knotY, knotZ;
        std::vector<Interp1D> derivs2;
};


#endif // end of include guard: INTERPOLATION_HH
