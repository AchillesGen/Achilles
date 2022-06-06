#ifndef PYBINDINGS_HH
#define PYBINDINGS_HH

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"
#include "pybind11/functional.h"

namespace py = pybind11;

template<typename... Args>
using overload_cast_ = py::detail::overload_cast_impl<Args...>;

// Utilities
void LoggingModule(py::module&);
void ConstantsModule(py::module&);
void InterpolationModule(py::module&);

// Physics Objects
void VectorModule(py::module&);
void ParticleInfoModule(py::module&);
void ParticleModule(py::module&);
void NucleusModule(py::module&);

// Calculation Objects
void InteractionsModule(py::module&);
void CascadeModule(py::module&);


#endif
