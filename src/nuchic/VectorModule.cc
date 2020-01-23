#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"

#include "nuchic/ThreeVector.hh"
#include "nuchic/FourVector.hh"

namespace py = pybind11;

PYBIND11_MODULE(vectors, m) {
    py::class_<ThreeVector, std::shared_ptr<ThreeVector>>(m, "Vector3")
        .def(py::init<>())
        .def(py::init<const std::array<double, 3>&>())
        .def(py::init<const double&, const double&, const double&>())
        .def(py::init<const ThreeVector&>())
        // Setters
        .def("set_pxpypz", (void (ThreeVector::*)(const std::array<double, 3>)) &ThreeVector::SetPxPyPz)
        .def("set_pxpypz", (void (ThreeVector::*)(const double&, const double&, const double&)) &ThreeVector::SetPxPyPz)
        .def("set_px", &ThreeVector::SetPx)
        .def("set_py", &ThreeVector::SetPy)
        .def("set_pz", &ThreeVector::SetPz)
        // Getters
        .def("position", &ThreeVector::Position)
        .def("px", &ThreeVector::Px)
        .def("py", &ThreeVector::Py)
        .def("pz", &ThreeVector::Pz)
        .def("pt2", &ThreeVector::Pt2)
        .def("pt", &ThreeVector::Pt)
        .def("p2", &ThreeVector::P2)
        .def("p", &ThreeVector::P)
        .def("magnitude2", &ThreeVector::Magnitude2)
        .def("magnitude", &ThreeVector::Magnitude)
        .def("theta", &ThreeVector::Theta)
        .def("phi", &ThreeVector::Phi)
        // Functions
        .def("dot", &ThreeVector::Dot)
        .def("cross", &ThreeVector::Cross)
        .def("unit", &ThreeVector::Unit)
        // Operator Overloads
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(-py::self)
        .def(py::self -= py::self)
        .def(float() * py::self)
        .def(py::self * float())
        .def(py::self *= float())
        .def(py::self / float())
        .def(py::self /= float())
        // Comparison Operators
        .def(py::self == py::self)
        .def(py::self != py::self)
        // Access Operators
        .def("__setitem__", (double& (ThreeVector::*)(const std::size_t&)) &ThreeVector::operator[])
        .def("__getitem__", (double& (ThreeVector::*)(const std::size_t&)) &ThreeVector::operator[])
        // String methods
        .def("__str__", &ThreeVector::ToString)
        .def("__repr__", &ThreeVector::ToString);

    py::class_<FourVector, std::shared_ptr<FourVector>>(m, "Vector4")
        .def(py::init<>())
        .def(py::init<const std::array<double, 4>&>())
        .def(py::init<const double&, const double&, const double&, const double&>())
        .def(py::init<const FourVector&>())
        // Setters
        .def("set_pxpypze", (void (FourVector::*)(const std::array<double, 4>)) &FourVector::SetPxPyPzE)
        .def("set_pxpypze", (void (FourVector::*)(const double&, const double&, const double&, const double&))
                &FourVector::SetPxPyPzE)
        .def("set_px", &FourVector::SetPx)
        .def("set_py", &FourVector::SetPy)
        .def("set_pz", &FourVector::SetPz)
        .def("set_e", &FourVector::SetE)
        // Getters
        .def("px", &FourVector::Px)
        .def("py", &FourVector::Py)
        .def("pz", &FourVector::Pz)
        .def("pt2", &FourVector::Pt2)
        .def("pt", &FourVector::Pt)
        .def("p2", &FourVector::P2)
        .def("p", &FourVector::P)
        .def("magnitude2", &FourVector::Magnitude2)
        .def("magnitude", &FourVector::Magnitude)
        .def("theta", &FourVector::Theta)
        .def("phi", &FourVector::Phi)
        .def("vec3", &FourVector::Vec3)
        // Functions
        .def("dot", &FourVector::Dot)
        .def("boost", (FourVector (FourVector::*)(const ThreeVector&)) &FourVector::Boost)
        .def("boost", (FourVector (FourVector::*)(const double&, const double&, const double&)) &FourVector::Boost)
        .def("boost_vector", &FourVector::BoostVector)
        // Operator Overloads
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(-py::self)
        .def(py::self -= py::self)
        .def(float() * py::self)
        .def(py::self * float())
        .def(py::self *= float())
        .def(py::self / float())
        .def(py::self /= float())
        // Comparison Operators
        .def(py::self == py::self)
        .def(py::self != py::self)
        // Access Operators
        .def("__setitem__", (double& (FourVector::*)(const std::size_t&)) &FourVector::operator[])
        .def("__getitem__", (double& (FourVector::*)(const std::size_t&)) &FourVector::operator[])
        // String methods
        .def("__str__", &FourVector::ToString)
        .def("__repr__", &FourVector::ToString);
}
