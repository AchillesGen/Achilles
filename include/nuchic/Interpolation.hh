#ifndef INTERPOLATION_HH
#define INTERPOLATION_HH

#include <vector>

// #include "pybind11/numpy.h"

// namespace py = pybind11;

// using pyArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

namespace nuchic {

enum class InterpolationType {
    NearestNeighbor,
    Polynomial,
    CubicSpline
};

double Polint(const std::vector<double>&, const std::vector<double>&,
              size_t, double);

/// Class to perform one-dimensional interpolations of data. Currently, only Cubic Splines are
/// implemented as an interpolator. The Cubic Spline is based off of the algorithm provided by
/// Numerical Recipes.
class Interp1D {
    public:
        /// @name Constructor and Destructor
        ///@{

        /// Constructor
        Interp1D() = default;
        Interp1D(const std::vector<double>&, const std::vector<double>&,
                 InterpolationType=InterpolationType::CubicSpline);
        Interp1D(const Interp1D&) = default;
        Interp1D(Interp1D&&) = default;
        Interp1D& operator=(const Interp1D&) = default;
        Interp1D& operator=(Interp1D&&) = default;

        /// Destructor
        ~Interp1D() = default;
        ///@}

        ///@name Interpolation Functions
        ///@{

        /// Create the cubic spline knots and derivatives for later interpolation
        ///@param x: A vector containing the x-components to interpolate over
        ///@param y: A vector containing the function values at the x-points
        ///@param derivLeft: The derivative of the left most knot
        ///@param derivRight: The derivative of the right most knot
        void CubicSpline(const double& derivLeft=maxDeriv, const double& derivRight=maxDeriv);

        const double& min() const { return knotX.front(); }
        const double& max() const { return knotX.back(); }

        void SetData(const std::vector<double> &x, const std::vector<double> &y) { knotX = x; knotY = y; }
        void SetType(InterpolationType mode) { kMode = mode; }
        void SetPolyOrder(size_t order) { polyOrder = order+1; }

        /// Function to perform the interpolation at the given input point
        ///@param x: Value to interpolate the function at
        ///@return double: The interpolated value of the function
        double operator()(const double&) const;
        ///@}

    private:
        double PolynomialInterp(double) const;

        InterpolationType kMode{InterpolationType::CubicSpline};
        static constexpr double maxDeriv = 1.E30;
        bool kSplineInit{};
        size_t polyOrder{4};
        std::vector<double> knotX, knotY, derivs2;
};

/// Class to perform two-dimensional interpolations of data. Currently, only Bicubic Splines are
/// implemented as an interpolator. The Bicubic Spline is based off of the algorithm provided by
/// Numerical Recipes.
class Interp2D {
    public:
        /// @name Constructor and Destructor
        ///@{

        /// Constructor
        Interp2D() = default;
        Interp2D(const std::vector<double>&, const std::vector<double>&,
                 const std::vector<double>&,
                 InterpolationType=InterpolationType::CubicSpline);
        // Interp2D(const std::vector<double>&, const std::vector<double>&,
        //          const pyArray&,
        //          InterpolationType=InterpolationType::CubicSpline);
        Interp2D(const Interp2D&) = default;
        Interp2D(Interp2D&&) = default;
        Interp2D& operator=(const Interp2D&) = default;
        Interp2D& operator=(Interp2D&&) = default;

        /// Destructor
        ~Interp2D() = default;
        ///@}

        ///@name Interpolation Functions
        ///@{

        /// Create the bicubic spline knots and derivatives for later interpolation
        ///@param x: A vector containing the x-components to interpolate over
        ///@param y: A vector containing the y-components to interpolate over
        ///@param z: A vector containing the function values at the x,y-points. The shape of
        ///          z should be given by x.size()*y.size() 
        void BicubicSpline();

        const double& Xmin() const { return knotX.front(); }
        const double& Xmax() const { return knotX.back(); }

        const double& Ymin() const { return knotY.front(); }
        const double& Ymax() const { return knotY.back(); }

        void SetData(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z) {
            knotX = x; knotY = y; knotZ = z;
        }
        void SetType(InterpolationType mode) { kMode = mode; }
        void SetPolyOrder(size_t orderX, size_t orderY) { 
            polyOrderX = orderX+1; polyOrderY = orderY+1;
        }

        /// Function to perform the interpolation at the given input point
        ///@param x: x-value to interpolate the function at
        ///@param y: y-value to interpolate the fucntion at
        ///@return double: The interpolated value of the function
        double operator()(const double&, const double&) const;
        ///@}

    private:
        double NearestNeighbor(double, double) const;
        double PolynomialInterp(double, double) const;

        bool kSplineInit{};
        InterpolationType kMode{InterpolationType::CubicSpline};
        size_t polyOrderX{4}, polyOrderY{4};
        std::vector<double> knotX, knotY, knotZ;
        std::vector<Interp1D> derivs2;
};

}

#endif // end of include guard: INTERPOLATION_HH
