#include "pybind11/pybind11.h"
#include "pybind11/operators.h"
#include "pybind11/stl.h"

#include "nuchic/ParticleInfo.hh"

namespace py = pybind11;

// These are for convenience
using nuchic::PID;
using nuchic::ParticleInfoEntry;
using nuchic::ParticleInfo;

PYBIND11_MODULE(particle_info, m) {
    py::class_<PID, std::shared_ptr<PID>> pid(m, "PID");
    py::class_<ParticleInfoEntry, std::shared_ptr<ParticleInfoEntry>> particleEntry(m, "ParticleInfoEntry");
    py::class_<ParticleInfo, std::shared_ptr<ParticleInfo>> particleInfo(m, "ParticleInfo");

    pid.def(py::init<const long int&>())
       .def("as_int", &PID::AsInt)
       .def("anti", &PID::Anti)
       .def(py::self == py::self)
       .def(py::self != py::self)
       .def(py::self < py::self)
       .def(py::self > py::self)
       // Undefined
       .def_property_readonly_static("undefined", &PID::undefined)
       // Quarks
       .def_property_readonly_static("down", [](py::object){ return PID::down(); })
       .def_property_readonly_static("up", [](py::object){ return PID::up(); })
       .def_property_readonly_static("strange", [](py::object){ return PID::strange(); })
       .def_property_readonly_static("charm", [](py::object){ return PID::charm(); })
       .def_property_readonly_static("bottom", [](py::object){ return PID::bottom(); })
       .def_property_readonly_static("top", [](py::object){ return PID::top(); })
       // Leptons
       .def_property_readonly_static("electron", [](py::object){ return PID::electron(); })
       .def_property_readonly_static("nu_electron", [](py::object){ return PID::nu_electron(); })
       .def_property_readonly_static("muon", [](py::object){ return PID::muon(); })
       .def_property_readonly_static("nu_muon",[](py::object){ return PID::nu_muon(); })
       .def_property_readonly_static("tau", [](py::object){ return PID::tau(); })
       .def_property_readonly_static("nu_tau", [](py::object){ return PID::nu_tau(); })
       // Gauge Bosons
       .def_property_readonly_static("gluon", [](py::object){ return PID::gluon(); })
       .def_property_readonly_static("photon", [](py::object){ return PID::photon(); })
       .def_property_readonly_static("Zboson", [](py::object){ return PID::Zboson(); })
       .def_property_readonly_static("Wboson", [](py::object){ return PID::Wboson(); })
       .def_property_readonly_static("Higgs", [](py::object){ return PID::Higgs(); })
       // Mesons
       .def_property_readonly_static("pion0", [](py::object){ return PID::pion0(); })
       .def_property_readonly_static("pionp", [](py::object){ return PID::pionp(); })
       // Baryons
       .def_property_readonly_static("proton", [](py::object){ return PID::proton(); })
       .def_property_readonly_static("neutron", [](py::object){ return PID::neutron(); })
       // Strings
       .def("__str__", [](const PID &id){ return fmt::format("{}", id.AsInt()); } )
       // Hash
       .def("__hash__", &PID::AsInt);

    particleEntry.def(py::init<>())
                 .def(py::init<const PID&, const double&, const double&, const int&, const int&,
                               const int&, const int&, const int&, const bool&, const bool&,
                               const std::string&, const std::string&>())
                 .def(py::init<const ParticleInfoEntry&>())
                 .def(py::self == py::self)
                 .def(py::self != py::self)
                 .def("__str__", &ParticleInfoEntry::ToString);

    particleInfo.def(py::init<const std::shared_ptr<ParticleInfoEntry>&, const bool&>())
                .def(py::init<const long int&>())
                .def(py::init<const PID&, const bool&>())
                .def(py::init<const ParticleInfo&>())
                .def("name", &ParticleInfo::Name)
                .def("id", &ParticleInfo::ID)
                .def("is_baryon", &ParticleInfo::IsBaryon)
                .def("is_hadron", &ParticleInfo::IsHadron)
                .def("is_bhadron", &ParticleInfo::IsBHadron)
                .def("is_chadron", &ParticleInfo::IsCHadron)
                .def("is_anti", &ParticleInfo::IsAnti)
                .def("is_fermion", &ParticleInfo::IsFermion)
                .def("is_boson", &ParticleInfo::IsBoson)
                .def("is_scalar", &ParticleInfo::IsScalar)
                .def("is_vector", &ParticleInfo::IsVector)
                .def("is_tensor", &ParticleInfo::IsTensor)
                .def("is_photon", &ParticleInfo::IsPhoton)
                .def("is_lepton", &ParticleInfo::IsLepton)
                .def("is_quark", &ParticleInfo::IsQuark)
                .def("is_gluon", &ParticleInfo::IsGluon)
                .def("int_charge", &ParticleInfo::IntCharge)
                .def("charge", &ParticleInfo::Charge)
                .def("int_spin", &ParticleInfo::IntSpin)
                .def("spin", &ParticleInfo::Spin)
                .def("self_anti", &ParticleInfo::SelfAnti)
                .def("majorana", &ParticleInfo::Majorana)
                .def("stable", &ParticleInfo::Stable)
                .def("is_stable", &ParticleInfo::IsStable)
                .def("is_massive", &ParticleInfo::IsMassive)
                .def("mass", &ParticleInfo::Mass)
                .def("width", &ParticleInfo::Width)
                .def("generate_lifetime", &ParticleInfo::GenerateLifeTime)
                .def(py::self == py::self)
                .def(py::self != py::self)
                .def_static("init_database", &ParticleInfo::InitDatabase)
                .def_static("print_database", &ParticleInfo::PrintDatabase)
                .def_static("database", &ParticleInfo::Database);
}
