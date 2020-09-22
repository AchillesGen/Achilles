#include "nuchic/PyBindings.hh"
#include "nuchic/ThreeVector.hh"
#include "nuchic/FourVector.hh"

// Convience names
using nuchic::ThreeVector;
using nuchic::FourVector;

void VectorModule(py::module &m) {
    py::class_<ThreeVector, std::shared_ptr<ThreeVector>>(m, "Vector3", py::module_local())
        .def(py::init<>())
        .def(py::init<const std::array<double, 3>&>())
        .def(py::init<const double&, const double&, const double&>())
        .def(py::init<const ThreeVector&>())
        // Setters
        .def("set_xyz", overload_cast_<const std::array<double, 3>&>()(&ThreeVector::SetXYZ))
        .def("set_xyz", overload_cast_<const double&, const double&, const double&>()(&ThreeVector::SetXYZ))
        .def("set_pxpypz", overload_cast_<const std::array<double, 3>>()(&ThreeVector::SetPxPyPz))
        .def("set_pxpypz", overload_cast_<const double&, const double&, const double&>()(&ThreeVector::SetPxPyPz))
        .def("set_x", &ThreeVector::SetX)
        .def("set_y", &ThreeVector::SetY)
        .def("set_z", &ThreeVector::SetZ)
        .def("set_px", &ThreeVector::SetPx)
        .def("set_py", &ThreeVector::SetPy)
        .def("set_pz", &ThreeVector::SetPz)
        // Getters
        .def("position", &ThreeVector::Position)
        .def("x", &ThreeVector::X)
        .def("y", &ThreeVector::Y)
        .def("z", &ThreeVector::Z)
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
        .def(py::self * py::self)
        .def(-py::self)
        // NOTE: There is a clang bug that creates a warning for the line below
        //.def(py::self -= py::self)
        // Therefore, we have to do the following as a workaround
        .def("__isub__", [](ThreeVector &a, const ThreeVector &b) {
                return a -= b;
            }, py::is_operator())
        .def(float() * py::self)
        .def(py::self * float())
        .def(py::self *= float())
        .def(py::self / float())
        .def(py::self /= float())
        // Comparison Operators
        .def(py::self == py::self)
        .def(py::self != py::self)
        // Access Operators
        .def("__setitem__", overload_cast_<const std::size_t&>()(&ThreeVector::operator[]))
        .def("__getitem__", overload_cast_<const std::size_t&>()(&ThreeVector::operator[]))
        // String methods
        .def("__str__", &ThreeVector::ToString)
        .def("__repr__", &ThreeVector::ToString);

    py::class_<FourVector, std::shared_ptr<FourVector>>(m, "Vector4", py::module_local())
        .def(py::init<>())
        .def(py::init<const std::array<double, 4>&>())
        .def(py::init<const double&, const double&, const double&, const double&>())
        .def(py::init<const ThreeVector&, const double&>())
        .def(py::init<const FourVector&>())
        // Setters
        .def("set_pxpypze", overload_cast_<const std::array<double, 4>>()(&FourVector::SetPxPyPzE))
        .def("set_pxpypze", overload_cast_<const double&, const double&, const double&, const double&>()(&FourVector::SetPxPyPzE))
        .def("set_px", &FourVector::SetPx)
        .def("set_py", &FourVector::SetPy)
        .def("set_pz", &FourVector::SetPz)
        .def("set_e", &FourVector::SetE)
        // Getters
        .def("x", &FourVector::X)
        .def("y", &FourVector::Y)
        .def("z", &FourVector::Z)
        .def("time", &FourVector::T)
        .def("px", &FourVector::Px)
        .def("py", &FourVector::Py)
        .def("pz", &FourVector::Pz)
        .def("energy", &FourVector::E)
        .def("pt2", &FourVector::Pt2)
        .def("pt", &FourVector::Pt)
        .def("p2", &FourVector::P2)
        .def("p", &FourVector::P)
        .def("mass", &FourVector::M)
        .def("mass2", &FourVector::M2)
        .def("magnitude2", &FourVector::Magnitude2)
        .def("magnitude", &FourVector::Magnitude)
        .def("theta", &FourVector::Theta)
        .def("phi", &FourVector::Phi)
        .def("vec3", &FourVector::Vec3)
        // Functions
        .def("dot", &FourVector::Dot)
        .def("boost", overload_cast_<const ThreeVector&>()(&FourVector::Boost))
        .def("boost", overload_cast_<const double&, const double&, const double&>()(&FourVector::Boost))
        .def("cross", &FourVector::Cross)
        .def("boost_vector", &FourVector::BoostVector)
        // Operator Overloads
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(-py::self)
        // NOTE: There is a clang bug that creates a warning for the line below
        //.def(py::self -= py::self)
        // Therefore, we have to do the following as a workaround
        .def("__isub__", [](FourVector &a, const FourVector &b) {
                return a -= b;
            }, py::is_operator())
        .def(float() * py::self)
        .def(py::self * float())
        .def(py::self *= float())
        .def(py::self / float())
        .def(py::self /= float())
        // Comparison Operators
        .def(py::self == py::self)
        .def(py::self != py::self)
        // Access Operators
        .def("__setitem__", overload_cast_<const std::size_t&>()(&FourVector::operator[]))
        .def("__getitem__", overload_cast_<const std::size_t&>()(&FourVector::operator[]))
        // String methods
        .def("__str__", &FourVector::ToString)
        .def("__repr__", &FourVector::ToString);
}
