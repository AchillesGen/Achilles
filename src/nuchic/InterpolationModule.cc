#include "nuchic/PyBindings.hh"
#include "nuchic/Interpolation.hh"

// These are for convenience
using nuchic::Interp1D;
using nuchic::Interp2D;

void Interp2D::BicubicSpline(const std::vector<double>& x, const std::vector<double>& y,
        const pyArray& z) {
    // Ensure inputs are valid
    if(!std::is_sorted(x.begin(), x.end()))
        std::runtime_error("Inputs must be increasing.");
    if(std::adjacent_find(x.begin(), x.end()) != x.end())
        std::runtime_error("Inputs must all be unique.");
    if(!std::is_sorted(y.begin(), y.end()))
        std::runtime_error("Inputs must be increasing.");
    if(std::adjacent_find(y.begin(), y.end()) != y.end())
        std::runtime_error("Inputs must all be unique.");

    // Ensure shape is correct
    if(z.ndim() != 2)
        std::runtime_error("Data should be in a 2D array");
    if(x.size() != static_cast<std::size_t>(z.shape()[0]))
        std::runtime_error("Input and output arrays must be the same size.");
    if(y.size() != static_cast<std::size_t>(z.shape()[1]))
        std::runtime_error("Input and output arrays must be the same size.");

    knotX = x;
    knotY = y;

    // Copy numpy array into std::vector<double>
    knotZ.resize(static_cast<std::size_t>(z.size()));
    std::memcpy(knotZ.data(), z.data(), static_cast<std::size_t>(z.size())*sizeof(double));

    for(std::size_t i = 0; i < x.size(); ++i) {
        derivs2.emplace_back();
        derivs2[i].CubicSpline(y, std::vector<double>(static_cast<std::size_t>(z.at(i, 0)), z.at(i, y.size()-1)));
    }

    kInit = true;
}

void InterpolationModule(py::module &m) {
    py::class_<Interp1D>(m, "Interp1D")
        // Constructors
        .def(py::init<>())
        // Functions
        .def("cubic_spline", &Interp1D::CubicSpline,
                py::arg("x"), py::arg("y"),
                py::arg("derivLeft") = 1e30, py::arg("derivRight") = 1e30)
        .def("call", &Interp1D::operator())
        .def("__call__", &Interp1D::operator());

    py::class_<Interp2D>(m, "Interp2D")
        // Constructors
        .def(py::init<>())
        // Functions
        .def("bicubic_spline",
                overload_cast_<const std::vector<double>&, const std::vector<double>&,
                               const std::vector<double>&>()(&Interp2D::BicubicSpline))
        .def("bicubic_spline",
                overload_cast_<const std::vector<double>&, const std::vector<double>&,
                               const pyArray&>()(&Interp2D::BicubicSpline))
        .def("call", &Interp2D::operator())
        .def("__call__", &Interp2D::operator());
}
