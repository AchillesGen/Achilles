#ifndef INTERPOLATION_HH
#define INTERPOLATION_HH

#include <vector>

#include "pybind11/numpy.h"

namespace py = pybind11;

using pyArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

namespace nuchic {

/// Class to perform one-dimensional interpolations of data. Currently, only Cubic Splines are
/// implemented as an interpolator. The Cubic Spline is based off of the algorithm provided by
/// Numerical Recipes.
class Interp1D {
    public:
        /// @name Constructor and Destructor
        ///@{

        /// Constructor
        Interp1D() : kInit(false) {}
        /// Destructor
        ~Interp1D() {}
        ///@}

        ///@name Interpolation Functions
        ///@{

        /// Create the cubic spline knots and derivatives for later interpolation
        ///@param x: A vector containing the x-components to interpolate over
        ///@param y: A vector containing the function values at the x-points
        ///@param derivLeft: The derivative of the left most knot
        ///@param derivRight: The derivative of the right most knot
        void CubicSpline(const std::vector<double>&, const std::vector<double>&,
                const double& derivLeft=1.e30, const double& derivRight=1.e30);

        /// Function to perform the interpolation at the given input point
        ///@param x: Value to interpolate the function at
        ///@return double: The interpolated value of the function
        double operator()(const double&) const;
        ///@}

    private:
        bool kInit;
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
        Interp2D() : kInit(false) {}
        /// Destructor
        ~Interp2D() {}
        ///@}

        ///@name Interpolation Functions
        ///@{

        /// Create the bicubic spline knots and derivatives for later interpolation
        ///@param x: A vector containing the x-components to interpolate over
        ///@param y: A vector containing the y-components to interpolate over
        ///@param z: A vector containing the function values at the x,y-points. The shape of
        ///          z should be given by x.size()*y.size() 
        void BicubicSpline(const std::vector<double>&, const std::vector<double>&,
                const std::vector<double>&);

        /// Create the bicubic spline knots and derivatives for later interpolation. This version
        /// is used to nicely interface with a 2D numpy array inside python.
        ///@param x: A vector containing the x-components to interpolate over
        ///@param y: A vector containing the y-components to interpolate over
        ///@param z: A numpy array containing the function values at the x,y-points. 
        ///          The shape of z should be given by (x.size(), y.size())
        void BicubicSpline(const std::vector<double>&, const std::vector<double>&,
                const pyArray&);

        /// Function to perform the interpolation at the given input point
        ///@param x: x-value to interpolate the function at
        ///@param y: y-value to interpolate the fucntion at
        ///@return double: The interpolated value of the function
        double operator()(const double&, const double&) const;
        ///@}

    private:
        bool kInit;
        std::vector<double> knotX, knotY, knotZ;
        std::vector<Interp1D> derivs2;
};

}

#endif // end of include guard: INTERPOLATION_HH
