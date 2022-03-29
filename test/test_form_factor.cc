#include "catch2/catch.hpp"
#include "catch2/trompeloeil.hpp"

#include <iostream>

#include "nuchic/FormFactor.hh"
#include "nuchic/Units.hh"

#include "nuchic/Utilities.hh"
#include "yaml-cpp/yaml.h"

using nuchic::operator""_GeV;

// class DummyFF : public nuchic::FormFactorImpl, nuchic::RegistrableFormFactor<DummyFF> {
//     public:
//         DummyFF(const YAML::Node&) {}
//         void Evaluate(double, nuchic::FormFactor::Values&) const override {}
//         nuchic::FFType Type() const override { return nuchic::FFType::vector; }
// 
//         static std::unique_ptr<nuchic::FormFactorImpl> Construct(const YAML::Node &node) {
//             return std::make_unique<DummyFF>(node);
//         }
//         static std::string Name() { return "MockFF"; }
// };

// class MockFF : public trompeloeil::mock_interface<DummyFF> {
//     static constexpr bool trompeloeil_movable_mock = true;
//     IMPLEMENT_CONST_MOCK2(Evaluate);
//     IMPLEMENT_CONST_MOCK0(Type);
// 
//     static std::unique_ptr<nuchic::FormFactorImpl> Construct(const YAML::Node&) {
//         return std::make_unique<MockFF>();
//     }
//     static std::string Name() { return "MockFF"; }
// };

void PrintFFVals(double q2, const nuchic::FormFactor::Values &vals) {
    std::cout << q2 << ",";
    std::cout << vals.F1p << "," << vals.F2p << ",";
    std::cout << vals.F1n << "," << vals.F2n << ",";
    std::cout << vals.FA << "," << vals.Fcoh << "\n";
}

void TestFF(nuchic::FormFactorImpl *ff) {
    static const double q2min = std::log10(10*10);
    static const double q2max = std::log10(5000*5000);
    static constexpr size_t nq2 = 1000;
    std::vector<double> q2_vec = nuchic::Logspace(q2min, q2max, nq2);
    for(const auto &q2 : q2_vec) {
        nuchic::FormFactor::Values vals;
        ff -> Evaluate(q2/1.0_GeV/1.0_GeV, vals);
        PrintFFVals(q2, vals);
    }
}

TEST_CASE("Factory", "[FormFactor]") {
}

TEST_CASE("VectorDipole", "[FormFactor]") {
    YAML::Node node = YAML::Load(R"node(
    lambda: 0.847
    Mu Proton: 2.79278
    Mu Neutron: -1.91315
    )node");
    nuchic::VectorDipole ff(node);
    std::cout << ff.Name() << "\n";
    TestFF(&ff);
}

TEST_CASE("KellyVector", "[FormFactor]") {
    YAML::Node node = YAML::Load(R"node(
    lambdasq: 0.7174
    Mu Proton: 2.79278
    Mu Neutron: -1.91315
    Gep Params: [-0.24, 10.98, 12.82, 21.97]
    Gen Params: [1.70, 3.30]
    Gmp Params: [0.12, 10.97, 18.86, 6.55]
    Gmn Params: [2.33, 14.72, 24.20, 84.1]
    )node");
    nuchic::Kelly ff(node);
    std::cout << ff.Name() << "\n";
    TestFF(&ff);
}

TEST_CASE("BBBAVector", "[FormFactor]") {
    YAML::Node node = YAML::Load(R"node(
    Mu Proton: 2.79278
    Mu Neutron: -1.91315
    NumeratorEp Params: [1.0, -0.0578, 0.0, 0.0]
    DenominatorEp Params: [11.1, 13.6, 33.0, 0.0]
    NumeratorEn Params: [0.0, 1.25, 1.3, 0.0]
    DenominatorEn Params: [9.86, 305.0, -758.0, 802.0]
    NumeratorMp Params: [1.0, 0.015, 0.0, 0.0]
    DenominatorMp Params: [11.1, 19.6, 7.54, 0.0]
    NumeratorMn Params: [1.0, 1.81, 0.0, 0.0]
    DenominatorMn Params: [14.1, 20.70, 68.7, 0.0]
    )node");
    nuchic::BBBA ff(node);
    std::cout << ff.Name() << "\n";
    TestFF(&ff);
}

TEST_CASE("ArringtonHillVector", "[FormFactor]") {
    YAML::Node node = YAML::Load(R"node(
    Mu Proton: 2.79278
    Mu Neutron: -1.91315
    tcut: 0.0779191396
    t0: -0.7
    Gep Params: [ 0.23916329807, -1.10985857441,  1.44438081306,  0.47956946560, 
                 -2.28689474187,  1.12663298498,  1.25061984354, -3.63102047159,
                  4.08221702379,  0.50409734650, -5.08512046051,  3.96774254395,
                 -0.98152907110 ]
    Gmp Params: [ 0.26414299414, -1.09530612212,  1.21855378178,  0.6611364935,
                 -1.40567892503, -1.35641843888,  1.44702915534,  4.2356697359,
                 -5.33404565341, -2.91630052096,  8.70740306757, -5.7069999438,
                  1.28081437589 ]
    Gen Params: [ 0.04891998138, -0.06452505391, -0.24082589738,  0.3921087449,
                  0.30044525860, -0.66188868718, -0.17563976969,  0.6246917245,
                 -0.07768429937, -0.23600397526,  0.09040197347,  0.0, 0.0 ]
    Gmn Params: [ 0.25775832696, -1.07954064206,  1.18218381220,  0.7110150858,
                 -1.34808093680, -1.66244402521,  2.62435442603,  1.7512344946,
                 -4.92230087889,  3.19789272731, -0.71207238995,  0.0, 0.0 ]
    )node");
    nuchic::ArringtonHill ff(node);
    std::cout << ff.Name() << "\n";
    TestFF(&ff);
}

TEST_CASE("AxialDipole", "[FormFactor]") {
    YAML::Node node = YAML::Load(R"node(
    MA: 1.000
    gan1: 1.2694
    gans: 0.08
    )node");
    nuchic::AxialDipole ff(node);
    std::cout << ff.Name() << "\n";
    TestFF(&ff);
}

TEST_CASE("HelmCoherent", "[FormFactor]") {
    YAML::Node node = YAML::Load(R"node(
    s: 1
    A: 12
    )node");
    nuchic::HelmFormFactor ff(node);
    std::cout << ff.Name() << "\n";
    TestFF(&ff);
}

TEST_CASE("LovatoCoherent", "[FormFactor]") {
    YAML::Node node = YAML::Load(R"node(
    b: 0.994885
    c: [6.01084, 0.143785, -4.06675, 1.45063, -0.135216]
    )node");
    nuchic::LovatoFormFactor ff(node);
    std::cout << ff.Name() << "\n";
    TestFF(&ff);
}
