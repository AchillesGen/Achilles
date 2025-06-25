#include "catch2/catch.hpp"

#include <iostream>

#include "Achilles/Constants.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/Units.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

using achilles::operator""_GeV;

class DummyFF : public achilles::FormFactorImpl, achilles::RegistrableFormFactor<DummyFF> {
  public:
    DummyFF() = default;
    void Evaluate(double, achilles::FormFactor::Values &vals) const override {
        vals.Gep = 1;
        vals.Gmp = 1;
        vals.Gen = 1;
        vals.Gmn = 1;
        Fill(1, vals);
    }

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(achilles::FFType, const YAML::Node &) {
        return std::make_unique<DummyFF>();
    }
    static std::string Name() { return "FFTest"; }
    static achilles::FFType Type() { return achilles::FFType::vector; }
};

TEST_CASE("Utilities", "[FormFactor]") {
    SECTION("Fill works") {
        DummyFF ff;
        achilles::FormFactor::Values vals;
        ff.Evaluate(1, vals);
        CHECK(vals.F1p == 1);
        CHECK(vals.F1n == 1);
        CHECK(vals.F2p == 0);
        CHECK(vals.F2n == 0);
    }
}

TEST_CASE("Vector", "[FormFactor]") {
    SECTION("Dipole") {
        YAML::Node node = YAML::Load("lambda: 1\nMu Proton: 1\nMu Neutron: 1");
        auto ff = achilles::VectorDipole::Construct(achilles::FFType::vector, node);
        CHECK_THROWS_WITH(achilles::VectorDipole::Construct(achilles::FFType::axial, node),
                          "FormFactor: Expected type vector, got type axial");
        achilles::FormFactor::Values vals;
        ff->Evaluate(1, vals);
        CHECK(vals.Gep == Approx(0.25));
        // Value below is the result of:
        // (-muN*Q2/(1+5.6*Q2/mp^2)/4mp^2) * Gep
        // with:
        // muN = 1, Q2 = 1, Gep = 1/4, mp = 0.93827208816
        CHECK(vals.Gen == Approx(-0.038578136359582302263831 * 0.25));
        CHECK(vals.Gmp == vals.Gep);
        CHECK(vals.Gmn == vals.Gep);
    }

    SECTION("Kelly") {
        YAML::Node node = YAML::Load(R"node(
            lambdasq: 1
            Mu Proton: 1
            Mu Neutron: 1
            Gep Params: [1, 1, 1, 1]
            Gen Params: [1, 1]
            Gmp Params: [1, 1, 1, 1]
            Gmn Params: [1, 1, 1, 1]
            )node");
        CHECK_THROWS_WITH(achilles::Kelly::Construct(achilles::FFType::axial, node),
                          "FormFactor: Expected type vector, got type axial");
        auto ff = achilles::Kelly::Construct(achilles::FFType::vector, node);
        achilles::FormFactor::Values vals;
        ff->Evaluate(1, vals);
        const double tau = 1.0 / (4 * pow(achilles::Constant::mp / 1_GeV, 2));
        CHECK(vals.Gep == Approx((1 + tau) / (1 + tau + tau * tau + tau * tau * tau)));
        CHECK(vals.Gen == Approx(0.25 * tau / (1 + tau)));
        CHECK(vals.Gmp == vals.Gep);
        CHECK(vals.Gmn == vals.Gep);
    }

    SECTION("BBBA") {
        YAML::Node node = YAML::Load(R"node(
            Mu Proton: 1
            Mu Neutron: 1
            NumeratorEp Params: [1, 1, 1, 1]
            DenominatorEp Params: [1, 1, 1, 1]
            NumeratorEn Params: [1, 1, 1, 1]
            DenominatorEn Params: [1, 1, 1, 1]
            NumeratorMp Params: [1, 1, 1, 1]
            DenominatorMp Params: [1, 1, 1, 1]
            NumeratorMn Params: [1, 1, 1, 1]
            DenominatorMn Params: [1, 1, 1, 1]
            )node");
        CHECK_THROWS_WITH(achilles::BBBA::Construct(achilles::FFType::axial, node),
                          "FormFactor: Expected type vector, got type axial");
        auto ff = achilles::BBBA::Construct(achilles::FFType::vector, node);
        achilles::FormFactor::Values vals;
        ff->Evaluate(1, vals);
        const double tau = 1.0 / (4 * pow(achilles::Constant::mp / 1_GeV, 2));
        const double num = 1 + tau + tau * tau + tau * tau * tau;
        const double den = 1 + tau + tau * tau + tau * tau * tau + tau * tau * tau * tau;
        CHECK(vals.Gep == Approx(num / den));
        CHECK(vals.Gen == Approx(num / den));
        CHECK(vals.Gmp == Approx(num / den));
        CHECK(vals.Gmn == Approx(num / den));
    }

    SECTION("Arrington-Hill") {
        YAML::Node node = YAML::Load(R"node(
            Mu Proton: 1
            Mu Neutron: 1
            tcut: 1
            t0: 1
            Gep Params: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            Gen Params: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            Gmp Params: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            Gmn Params: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            )node");
        CHECK_THROWS_WITH(achilles::ArringtonHill::Construct(achilles::FFType::axial, node),
                          "FormFactor: Expected type vector, got type axial");
        auto ff = achilles::ArringtonHill::Construct(achilles::FFType::vector, node);
        achilles::FormFactor::Values vals;
        ff->Evaluate(1, vals);
        CHECK(vals.Gep == Approx(13));
        CHECK(vals.Gen == Approx(13));
        CHECK(vals.Gmp == Approx(13));
        CHECK(vals.Gmn == Approx(13));
    }
}

TEST_CASE("Axial", "[FormFactor]") {
    SECTION("Dipole") {
        YAML::Node node = YAML::Load("MA: 1\ngan1: 1\ngans: 1");
        CHECK_THROWS_WITH(achilles::AxialDipole::Construct(achilles::FFType::vector, node),
                          "FormFactor: Expected type axial, got type vector");
        auto ff = achilles::AxialDipole::Construct(achilles::FFType::axial, node);
        achilles::FormFactor::Values vals;
        ff->Evaluate(1, vals);
        CHECK(vals.FA == Approx(-0.25));
        CHECK(vals.FAs == Approx(-0.25));
    }
    SECTION("AxialZExpansion") {
        YAML::Node node = YAML::Load(R"node(
            tcut: 0.1753180641  # [GeV^2] 9*mpi^2 = 9*(0.13957)^2
            t0: -0.28 # [GeV^2]
            CC Params: [-0.759, 2.30, -0.6, -3.8, 2.3, 2.16, -0.896, -1.58, 0.823]
            Strange Params: [1, 1, 1, 1, 1, 1, 1, 1, 1]
            )node");
        CHECK_THROWS_WITH(achilles::AxialZExpansion::Construct(achilles::FFType::vector, node),
                          "FormFactor: Expected type axial, got type vector");
        auto ff = achilles::AxialZExpansion::Construct(achilles::FFType::axial, node);
        achilles::FormFactor::Values vals;
        ff->Evaluate(0, vals);
        // fA(Q2=0) = gA = -1.27...
        CHECK(vals.FA == Approx(-1.27).epsilon(1e-2));
        // Q2=0 --> z=-0.234172
        // \sum_{n=0}^{8} (-0.234172)^n \approx 0.810262
        CHECK(vals.FAs == Approx(0.810262).epsilon(1e-6));
    }
}

TEST_CASE("Coherent", "[FormFactor]") {
    SECTION("Helm") {
        const double s = 0.2;
        const double A = 12;
        YAML::Node node = YAML::Load(fmt::format("s: {}\nA: {}", s, A));
        auto ff = achilles::HelmFormFactor::Construct(achilles::FFType::coherent, node);
        CHECK_THROWS_WITH(achilles::HelmFormFactor::Construct(achilles::FFType::vector, node),
                          "FormFactor: Expected type coherent, got type vector");
        achilles::FormFactor::Values vals;
        ff->Evaluate(1, vals);
        // Dummy result for A=12, s=0.2
        const double r = sqrt(1.2 * 1.2 * std::cbrt(A) * std::cbrt(A) - 5 * s * s);
        const double kappa = 1.0 / achilles::Constant::HBARC;
        const double result = 3 * exp(-kappa * kappa * 0.2 * 0.2 / 2) *
                              (sin(kappa * r) - kappa * r * cos(kappa * r)) / pow(kappa * r, 3);
        CHECK(vals.Fcoh == Approx(result));
    }

    SECTION("Lovato") {
        YAML::Node node = YAML::Load("b: 1\nc: [1, 1, 1, 1, 1]");
        CHECK_THROWS_WITH(achilles::LovatoFormFactor::Construct(achilles::FFType::vector, node),
                          "FormFactor: Expected type coherent, got type vector");
        auto ff = achilles::LovatoFormFactor::Construct(achilles::FFType::coherent, node);
        achilles::FormFactor::Values vals;
        ff->Evaluate(1, vals);
        const double x = 1.0 / achilles::Constant::HBARC;
        const double result = exp(-0.5 * x * x) * (1 + x + x * x + x * x * x + x * x * x * x) / 6.0;
        CHECK(vals.Fcoh == Approx(result));
    }
}

TEST_CASE("Builder", "[FormFactor]") {
    YAML::Node vector = YAML::Load("lambda: 1\nMu Proton: 1\nMu Neutron: 1");
    YAML::Node axial = YAML::Load("MA: 1\ngan1: 1\ngans: 1");
    YAML::Node coherent = YAML::Load("b: 1\nc: [1, 1, 1, 1, 1]");
    auto ff = achilles::FormFactorBuilder()
                  .Vector("VectorDipole", vector)
                  .AxialVector("AxialDipole", axial)
                  .Coherent("Lovato", coherent)
                  .build();
    auto vals = ff->operator()(1);
    const double x = 1.0 / achilles::Constant::HBARC;
    const double result = exp(-0.5 * x * x) * (1 + x + x * x + x * x * x + x * x * x * x) / 6.0;
    CHECK(vals.Gep == Approx(0.25));
    // Value below is the result of:
    // (-muN*Q2/(1+5.6*Q2/mp^2)/4mp^2) * Gep
    // with:
    // muN = 1, Q2 = 1, Gep = 1/4, mp = 0.93827208816
    CHECK(vals.Gen == Approx(-0.038578136359582302263831 * 0.25));
    CHECK(vals.Gmp == vals.Gep);
    CHECK(vals.Gmn == vals.Gep);
    CHECK(vals.FA == Approx(-0.25));
    CHECK(vals.FAs == Approx(-0.25));
    CHECK(vals.Fcoh == Approx(result));
}
