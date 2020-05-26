#include "pybind11/pybind11.h"

#include "nuchic/Constants.hh"

namespace py = pybind11;
namespace Constant = nuchic::Constant;

PYBIND11_MODULE(constants, m) {
    m.attr("c") = py::float_(Constant::C);
    m.attr("hbarc") = py::float_(Constant::HBARC);
    m.attr("hbarc2") = py::float_(Constant::HBARC2);
    m.attr("mp") = py::float_(Constant::mp);
    m.attr("mn") = py::float_(Constant::mn);
    m.attr("mN") = py::float_(Constant::mN);
}
