// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"

#include "Achilles/FormFactor.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Variations.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include <cmath>

using achilles::FormFactorVariation;

namespace {
const YAML::Node BaseConfig() {
    return YAML::Load(R"(
vector: VectorDipole
axial: AxialDipole
coherent: Helm
resonancevector: ResonanceVectorDummy
resonanceaxial: ResonanceAxialDummy
mecvector: MesonExchangeVector
mecaxial: MesonExchangeAxial
hyperon: Hyperon
VectorDipole: {lambda: 0.84, Mu Proton: 2.79278, Mu Neutron: -1.91315}
AxialDipole: {MA: 1.0, gan1: 1.2694, gans: 0.08}
Helm: {s: 1, A: 12}
ResonanceVectorDummy: {resV: 1}
ResonanceAxialDummy: {resA: 1}
MesonExchangeVector: {MvSq: 0.71, cv3norm: 2.13, cv4norm: -1.15, cv5norm: 0.48}
MesonExchangeAxial: {MaDeltaSq: 1.1025, ca5norm: 1.2}
Hyperon: {dummy: 0}
AxialZExpansion:
  tcut: 0.1753180641
  t0: -0.28
  CC Params: [-0.759, 2.30, -0.6, -3.8, 2.3, 2.16, -0.896, -1.58, 0.823]
  Strange Params: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
BBBA:
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
)");
}
} // namespace

TEST_CASE("Form factor config edits", "[Variations]") {
    auto base = BaseConfig();

    SECTION("Overrides are merged without modifying the nominal") {
        auto merged =
            FormFactorVariation::Merge(base, YAML::Load("{AxialDipole: {MA: 1.2}, vector: BBBA}"));
        CHECK(merged["AxialDipole"]["MA"].as<double>() == 1.2);
        CHECK(merged["AxialDipole"]["gan1"].as<double>() == 1.2694);
        CHECK(merged["vector"].as<std::string>() == "BBBA");
        CHECK(base["AxialDipole"]["MA"].as<double>() == 1.0);
        CHECK(base["vector"].as<std::string>() == "VectorDipole");
    }

    SECTION("Parameters are addressed by path, including sequence elements") {
        CHECK(FormFactorVariation::Lookup(base, "AxialDipole/MA").as<double>() == 1.0);
        CHECK(FormFactorVariation::Lookup(base, "AxialZExpansion/CC Params/1").as<double>() ==
              2.30);
        auto edited = FormFactorVariation::SetParameter(base, "AxialZExpansion/CC Params/1", 2.5);
        CHECK(FormFactorVariation::Lookup(edited, "AxialZExpansion/CC Params/1").as<double>() ==
              2.5);
        CHECK(FormFactorVariation::Lookup(edited, "AxialZExpansion/CC Params/2").as<double>() ==
              -0.6);
        CHECK(FormFactorVariation::Lookup(base, "AxialZExpansion/CC Params/1").as<double>() ==
              2.30);
        // Looking up a parameter must not insert it
        CHECK_THROWS(FormFactorVariation::Lookup(base, "AxialDipole/MV"));
        CHECK_FALSE(base["AxialDipole"]["MV"]);
        CHECK_THROWS(FormFactorVariation::Lookup(base, "AxialZExpansion/CC Params/20"));
        CHECK_THROWS(FormFactorVariation::SetParameter(base, "AxialZExpansion/CC Params", 1));
    }
}

TEST_CASE("Alternate form factors can be built from edited configs", "[Variations]") {
    auto base = BaseConfig();
    achilles::FormFactorBuilder builder;
    auto nominal = achilles::NuclearModel::BuildFormFactor(base, builder);

    SECTION("Parameter variation") {
        auto varied = achilles::NuclearModel::BuildFormFactor(
            FormFactorVariation::SetParameter(base, "AxialDipole/MA", 1.2), builder);
        const double Q2 = 0.5;
        const double expected = std::pow((1 + Q2 / 1.0) / (1 + Q2 / 1.44), 2);
        CHECK((*varied)(Q2).FA / (*nominal)(Q2).FA == Catch::Approx(expected));
        CHECK((*varied)(Q2).F1p == (*nominal)(Q2).F1p);
    }

    SECTION("Functional form variation") {
        auto varied = achilles::NuclearModel::BuildFormFactor(
            FormFactorVariation::Merge(base, YAML::Load("{vector: BBBA}")), builder);
        CHECK((*varied)(1.0).F1p != Catch::Approx((*nominal)(1.0).F1p));
        CHECK((*varied)(1.0).FA == (*nominal)(1.0).FA);
    }
}

TEST_CASE("Variation groups define central values and members", "[Variations]") {
    SECTION("MinMax") {
        achilles::FormFactorGroup group(YAML::Load(
            "{Name: MA, Type: FormFactor, Parameter: AxialDipole/MA, MinMax: [0.9, 1.1]}"));
        CHECK(group.WeightNames() == std::vector<std::string>{"MA:central", "MA:min", "MA:max"});
        CHECK(group.GetCombination().Type() == "Envelope");
        auto info = group.Info();
        CHECK(info.parameter == "AxialDipole/MA");
        CHECK(info.member_values == std::vector<double>{0.9, 1.1});
    }

    SECTION("Scan") {
        achilles::FormFactorGroup group(
            YAML::Load("{Name: MA, Type: FormFactor, Parameter: AxialDipole/MA, Central: 1.05,"
                       " Scan: {Min: 0.9, Max: 1.2, Steps: 4}}"));
        CHECK(group.WeightNames() ==
              std::vector<std::string>{"MA:central", "MA:0.9", "MA:1", "MA:1.1", "MA:1.2"});
        CHECK(group.GetCombination().Type() == "None");
    }

    SECTION("Values with an explicit combination") {
        achilles::FormFactorGroup group(YAML::Load(
            "{Name: MA, Type: FormFactor, Parameter: AxialDipole/MA, Values: [0.95, 1.05],"
            " Combination: {Type: SymmetricHessian, Scale: 0.5}}"));
        CHECK(group.GetCombination().Type() == "SymmetricHessian");
        CHECK(group.GetCombination().Scale() == 0.5);
    }

    SECTION("Alternatives") {
        achilles::FormFactorGroup group(YAML::Load(R"(
Name: VectorFF
Type: FormFactor
Alternatives:
  - {Name: BBBA, Overrides: {vector: BBBA}}
  - {Name: Lambda, Parameter: VectorDipole/lambda, Value: 0.9}
)"));
        CHECK(group.WeightNames() ==
              std::vector<std::string>{"VectorFF:central", "VectorFF:BBBA", "VectorFF:Lambda"});
        CHECK(group.GetCombination().Type() == "Envelope");
    }

    SECTION("Invalid definitions") {
        auto build = [](const std::string &yaml) {
            return achilles::FormFactorGroup(YAML::Load(yaml));
        };
        // No members
        CHECK_THROWS(build("{Name: MA, Type: FormFactor, Parameter: AxialDipole/MA}"));
        // More than one member specification
        CHECK_THROWS(build("{Name: MA, Type: FormFactor, Parameter: AxialDipole/MA,"
                           " MinMax: [0.9, 1.1], Values: [1.2]}"));
        CHECK_THROWS(build("{Name: MA, Type: FormFactor, Parameter: AxialDipole/MA,"
                           " MinMax: [1.1, 0.9]}"));
        CHECK_THROWS(build("{Name: M A, Type: FormFactor, Parameter: AxialDipole/MA,"
                           " MinMax: [0.9, 1.1]}"));
        CHECK_THROWS(build("{Name: MA, Type: FormFactor, Parameter: AxialDipole/MA,"
                           " Values: [0.9, 0.9]}"));
        // Hessian needs pairs
        CHECK_THROWS(build("{Name: MA, Type: FormFactor, Parameter: AxialDipole/MA,"
                           " Values: [0.9, 1.0, 1.1], Combination: Hessian}"));
        CHECK_THROWS(build("{Name: MA, Type: FormFactor, Parameter: AxialDipole/MA,"
                           " Values: [0.9], Combination: NotACombination}"));
    }
}

TEST_CASE("Combinations into uncertainty bands", "[Variations]") {
    const double central = 10;
    auto build = [](const std::string &type) {
        return achilles::Combination::Build(YAML::Load(type), "None");
    };

    SECTION("Envelope") {
        auto band = build("Envelope")->Combine(central, {9, 12, 10.5});
        CHECK(band.lower == 9);
        CHECK(band.upper == 12);
        // The central value is always inside the envelope
        band = build("Envelope")->Combine(central, {11, 12});
        CHECK(band.lower == 10);
    }

    SECTION("Symmetric Hessian") {
        auto band = build("SymmetricHessian")->Combine(central, {13, 6});
        CHECK(band.lower == Catch::Approx(5));
        CHECK(band.upper == Catch::Approx(15));
    }

    SECTION("Asymmetric Hessian") {
        // Eigenvector 1: (+3, -1), eigenvector 2: (+0, -4) both down
        auto band = build("Hessian")->Combine(central, {13, 9, 6, 8});
        CHECK(band.upper == Catch::Approx(13));
        CHECK(band.lower == Catch::Approx(10 - std::sqrt(1 + 16)));
    }

    SECTION("Replicas") {
        auto band = build("Replicas")->Combine(central, {9, 11, 9, 11});
        const double sigma = std::sqrt(4.0 / 3.0);
        CHECK(band.lower == Catch::Approx(10 - sigma));
        CHECK(band.upper == Catch::Approx(10 + sigma));
    }

    SECTION("Scale") {
        auto band = build("{Type: SymmetricHessian, Scale: 0.5}")->Combine(central, {14});
        CHECK(band.lower == Catch::Approx(8));
        CHECK(band.upper == Catch::Approx(12));
    }

    SECTION("None") {
        auto band = build("None")->Combine(central, {5, 15});
        CHECK(band.lower == central);
        CHECK(band.upper == central);
    }
}

TEST_CASE("VariationHandler builds groups through the factory", "[Variations]") {
    achilles::VariationHandler handler(YAML::Load(R"(
- {Name: MA, Type: FormFactor, Parameter: AxialDipole/MA, MinMax: [0.9, 1.1]}
- Name: VectorFF
  Type: FormFactor
  Alternatives: [{Name: BBBA, Overrides: {vector: BBBA}}]
)"));
    CHECK(handler.NWeights() == 5);
    CHECK(handler.WeightNames() == std::vector<std::string>{"MA:central", "MA:min", "MA:max",
                                                            "VectorFF:central", "VectorFF:BBBA"});

    CHECK_THROWS(achilles::VariationHandler(YAML::Load("[{Name: x, Type: NotAVariation}]")));
    CHECK_THROWS(achilles::VariationHandler(YAML::Load(R"(
- {Name: MA, Type: FormFactor, Parameter: AxialDipole/MA, MinMax: [0.9, 1.1]}
- {Name: MA, Type: FormFactor, Parameter: AxialDipole/MA, MinMax: [0.8, 1.2]}
)")));
}

TEST_CASE("Edits never modify the nominal configuration", "[Variations]") {
    const auto nominal = BaseConfig();
    const auto reference = YAML::Dump(nominal);

    FormFactorVariation::Edit parameter;
    parameter.parameters = {{"AxialDipole/MA", 0.9, false}};
    auto edited = FormFactorVariation::EditedConfig(nominal, parameter);
    CHECK(edited["AxialDipole"]["MA"].as<double>() == 0.9);

    FormFactorVariation::Edit overrides;
    overrides.overrides = YAML::Load("{vector: BBBA, AxialDipole: {gan1: 1.3}}");
    edited = FormFactorVariation::EditedConfig(nominal, overrides);
    CHECK(edited["vector"].as<std::string>() == "BBBA");
    CHECK(edited["AxialDipole"]["MA"].as<double>() == 1.0);

    FormFactorVariation::Edit both = overrides;
    both.parameters = {{"AxialDipole/MA", 1.2, false}};
    edited = FormFactorVariation::EditedConfig(nominal, both);
    CHECK(edited["AxialDipole"]["MA"].as<double>() == 1.2);
    CHECK(edited["AxialDipole"]["gan1"].as<double>() == 1.3);

    // Each edit starts from the untouched nominal
    CHECK(YAML::Dump(nominal) == reference);

    // A parameter set to its generation value is recognised as no change
    FormFactorVariation::Edit same;
    same.parameters = {{"AxialDipole/MA", 1.0, false}};
    CHECK(FormFactorVariation::EditedConfig(nominal, same).IsNull());
}

namespace {
void CheckEquivalent(const std::string &shorthand, const std::string &verbose) {
    achilles::VariationHandler compact(YAML::Load(shorthand));
    achilles::VariationHandler expanded(YAML::Load(verbose));
    CHECK(compact.WeightNames() == expanded.WeightNames());
    const auto info_compact = compact.Info();
    const auto info_expanded = expanded.Info();
    REQUIRE(info_compact.size() == info_expanded.size());
    for(size_t i = 0; i < info_compact.size(); ++i) {
        CHECK(info_compact[i].name == info_expanded[i].name);
        CHECK(info_compact[i].type == info_expanded[i].type);
        CHECK(info_compact[i].combination == info_expanded[i].combination);
        CHECK(info_compact[i].combination_scale == info_expanded[i].combination_scale);
        CHECK(info_compact[i].parameter == info_expanded[i].parameter);
        CHECK(info_compact[i].member_values == info_expanded[i].member_values);
        CHECK(info_compact[i].weight_names == info_expanded[i].weight_names);
    }
}
} // namespace

TEST_CASE("Shorthand variations match the explicit groups", "[Variations]") {
    SECTION("Central with up and down") {
        CheckEquivalent("{FormFactors: {AxialDipole/MA: 1.0 +0.3 -0.1}}", R"(
- {Name: AxialDipole.MA, Type: FormFactor, Parameter: AxialDipole/MA, Central: 1.0,
   MinMax: [0.9, 1.3]}
)");
        CheckEquivalent("{FormFactors: {AxialDipole/MA: 1.0 -0.1 +0.3}}", R"(
- {Name: AxialDipole.MA, Type: FormFactor, Parameter: AxialDipole/MA, Central: 1.0,
   MinMax: [0.9, 1.3]}
)");
        CheckEquivalent("{FormFactors: {AxialDipole/MA: 1.05 +- 0.1}}", R"(
- {Name: AxialDipole.MA, Type: FormFactor, Parameter: AxialDipole/MA, Central: 1.05,
   MinMax: [0.95, 1.15]}
)");
    }

    SECTION("Scans") {
        CheckEquivalent("{FormFactors: {Kelly/lambdasq: {Min: 0.6, Max: 0.8, Number: 5}}}", R"(
- {Name: Kelly.lambdasq, Type: FormFactor, Parameter: Kelly/lambdasq,
   Scan: {Min: 0.6, Max: 0.8, Steps: 5}}
)");
        CheckEquivalent("{FormFactors: {Kelly/lambdasq: {Min: 0.6, Max: 0.8, Step: 0.05}}}", R"(
- {Name: Kelly.lambdasq, Type: FormFactor, Parameter: Kelly/lambdasq,
   Scan: {Min: 0.6, Max: 0.8, Steps: 5}}
)");
        CheckEquivalent("{FormFactors: {AxialZExpansion/CC Params/1: [2.2, 2.4]}}", R"(
- {Name: AxialZExpansion.CC_Params.1, Type: FormFactor,
   Parameter: AxialZExpansion/CC Params/1, Values: [2.2, 2.4]}
)");
        CheckEquivalent("{FormFactors: {AxialDipole/MA: 1.2}}", R"(
- {Name: AxialDipole.MA, Type: FormFactor, Parameter: AxialDipole/MA, Values: [1.2]}
)");
    }

    SECTION("Explicit options inside the shorthand") {
        CheckEquivalent(
            "{FormFactors: {AxialDipole/MA: {Values: [1.1, 0.9], Combination: Hessian}}}", R"(
- {Name: AxialDipole.MA, Type: FormFactor, Parameter: AxialDipole/MA, Values: [1.1, 0.9],
   Combination: Hessian}
)");
    }

    SECTION("Functional forms and spectral functions") {
        CheckEquivalent(R"(
FormFactors: {vector: [BBBA, VectorDipole], axial: AxialZExpansion}
SpectralFunctions: {MF: [pke12p_MF.data, pke12n_MF.data]}
)",
                        R"(
- Name: vector
  Type: FormFactor
  Alternatives:
    - {Name: BBBA, Overrides: {vector: BBBA}}
    - {Name: VectorDipole, Overrides: {vector: VectorDipole}}
- Name: axial
  Type: FormFactor
  Alternatives: [{Name: AxialZExpansion, Overrides: {axial: AxialZExpansion}}]
- Name: SF
  Type: SpectralFunction
  Alternatives:
    - {Name: MF, SpectralP: data/Spectral_Functions/pke12p_MF.data,
       SpectralN: data/Spectral_Functions/pke12n_MF.data}
)");
    }

    SECTION("Weight names") {
        achilles::VariationHandler handler(YAML::Load(R"(
FormFactors:
  AxialDipole/MA: 1.0 +- 0.1
  vector: [BBBA]
)"));
        CHECK(handler.WeightNames() ==
              std::vector<std::string>{"AxialDipole.MA:central", "AxialDipole.MA:min",
                                       "AxialDipole.MA:max", "vector:central", "vector:BBBA"});
    }

    SECTION("Invalid shorthand") {
        auto build = [](const std::string &yaml) {
            return achilles::VariationHandler(YAML::Load(yaml));
        };
        CHECK_THROWS(build("{FromFactors: {AxialDipole/MA: [1.1]}}"));
        CHECK_THROWS(build("{FormFactors: {MA: [1.1]}}"));
        CHECK_THROWS(build("{FormFactors: {AxialDipole/MA: about one}}"));
        CHECK_THROWS(build("{FormFactors: {AxialDipole/MA: {Min: 0.9, Max: 1.1}}}"));
        CHECK_THROWS(build("{FormFactors: {AxialDipole/MA: {Min: 0.9, Max: 1.1, Step: 0.15}}}"));
        CHECK_THROWS(build("{FormFactors: {AxialDipole/MA: {Values: [1.1], Typo: 1}}}"));
        CHECK_THROWS(build("{SpectralFunctions: {MF: [pke12p_MF.data]}}"));
    }
}

namespace {
achilles::Matrix MultiplyTransposed(const achilles::Matrix &a) {
    // a * a^T
    const size_t n = a.size();
    achilles::Matrix result(n, std::vector<double>(n, 0));
    for(size_t i = 0; i < n; ++i)
        for(size_t j = 0; j < n; ++j)
            for(size_t k = 0; k < a[i].size(); ++k) result[i][j] += a[i][k] * a[j][k];
    return result;
}

achilles::Matrix FromShifts(const achilles::EigenShifts &shifts, size_t n) {
    // sum_k sigma_k^2 v_k v_k^T
    achilles::Matrix result(n, std::vector<double>(n, 0));
    for(size_t k = 0; k < shifts.sigmas.size(); ++k)
        for(size_t i = 0; i < n; ++i)
            for(size_t j = 0; j < n; ++j)
                result[i][j] += shifts.sigmas[k] * shifts.sigmas[k] * shifts.directions[k][i] *
                                shifts.directions[k][j];
    return result;
}

void CheckMatrix(const achilles::Matrix &a, const achilles::Matrix &b) {
    REQUIRE(a.size() == b.size());
    for(size_t i = 0; i < a.size(); ++i)
        for(size_t j = 0; j < a.size(); ++j) CHECK(a[i][j] == Catch::Approx(b[i][j]).margin(1e-12));
}

// A generic symmetric positive definite 5x5 matrix: L L^T with a lower triangular L
const achilles::Matrix kCovariance = MultiplyTransposed({{1.0, 0, 0, 0, 0},
                                                         {0.4, 0.8, 0, 0, 0},
                                                         {-0.3, 0.2, 0.5, 0, 0},
                                                         {0.1, -0.6, 0.3, 0.7, 0},
                                                         {0.2, 0.1, -0.4, 0.2, 0.3}});
} // namespace

TEST_CASE("Symmetric eigen decomposition", "[Variations]") {
    SECTION("2x2") {
        auto eigen = achilles::SymmetricEigen({{2, 1}, {1, 2}});
        CHECK(eigen.values[0] == Catch::Approx(3));
        CHECK(eigen.values[1] == Catch::Approx(1));
        CHECK(eigen.vectors[0][0] == Catch::Approx(std::sqrt(0.5)));
        CHECK(eigen.vectors[0][1] == Catch::Approx(std::sqrt(0.5)));
        CHECK(eigen.vectors[1][0] == Catch::Approx(std::sqrt(0.5)));
        CHECK(eigen.vectors[1][1] == Catch::Approx(-std::sqrt(0.5)));
    }

    SECTION("Reconstruction and orthonormality") {
        auto eigen = achilles::SymmetricEigen(kCovariance);
        const size_t n = kCovariance.size();
        achilles::Matrix rebuilt(n, std::vector<double>(n, 0));
        for(size_t k = 0; k < n; ++k) {
            if(k > 0) CHECK(eigen.values[k] <= eigen.values[k - 1]);
            for(size_t l = 0; l < n; ++l) {
                double dot = 0;
                for(size_t i = 0; i < n; ++i) dot += eigen.vectors[k][i] * eigen.vectors[l][i];
                CHECK(dot == Catch::Approx(k == l ? 1.0 : 0.0).margin(1e-12));
            }
            for(size_t i = 0; i < n; ++i)
                for(size_t j = 0; j < n; ++j)
                    rebuilt[i][j] += eigen.values[k] * eigen.vectors[k][i] * eigen.vectors[k][j];
        }
        CheckMatrix(rebuilt, kCovariance);
    }
}

TEST_CASE("Covariance inputs", "[Variations]") {
    SECTION("Correlation and uncertainties") {
        auto cov = achilles::ParseCovariance(
            YAML::Load("{Uncertainties: [0.1, 0.2], Correlation: [[1, -0.5], [-0.5, 1]]}"), 2);
        CheckMatrix(cov, {{0.01, -0.01}, {-0.01, 0.04}});
        cov = achilles::ParseCovariance(YAML::Load("{Uncertainties: [0.1, 0.2]}"), 2);
        CheckMatrix(cov, {{0.01, 0}, {0, 0.04}});
    }

    SECTION("Hessian is the inverse covariance") {
        YAML::Node node;
        node["Covariance"] = kCovariance;
        const auto cov = achilles::ParseCovariance(node, 5);
        // Invert through the eigen decomposition to build the Hessian input
        const auto eigen = achilles::SymmetricEigen(kCovariance);
        achilles::Matrix hessian(5, std::vector<double>(5, 0));
        for(size_t k = 0; k < 5; ++k)
            for(size_t i = 0; i < 5; ++i)
                for(size_t j = 0; j < 5; ++j)
                    hessian[i][j] += eigen.vectors[k][i] * eigen.vectors[k][j] / eigen.values[k];
        YAML::Node hess_node;
        hess_node["Hessian"] = hessian;
        CheckMatrix(achilles::ParseCovariance(hess_node, 5), cov);
    }

    SECTION("Invalid inputs") {
        auto parse = [](const std::string &yaml, size_t n) {
            return achilles::ParseCovariance(YAML::Load(yaml), n);
        };
        CHECK_THROWS(parse("{Uncertainties: [0.1, 0.2], Covariance: [[1, 0], [0, 1]]}", 2));
        CHECK_THROWS(parse("{Covariance: [[1, 0.5], [0.4, 1]]}", 2));
        CHECK_THROWS(parse("{Covariance: [[1, 0], [0, 1]]}", 3));
        CHECK_THROWS(parse("{Uncertainties: [0.1, 0.2], Correlation: [[1, 1.5], [1.5, 1]]}", 2));
        CHECK_THROWS(parse("{Uncertainties: [0.1, 0.2], Correlation: [[2, 0], [0, 1]]}", 2));
        CHECK_THROWS(parse("{Correlation: [[1, 0], [0, 1]]}", 2));
        CHECK_THROWS(parse("{Hessian: [[1, 0], [0, -1]]}", 2));
        CHECK_THROWS(parse("{}", 2));
    }
}

TEST_CASE("Eigenvector shifts", "[Variations]") {
    SECTION("Shifts reproduce the covariance") {
        auto shifts = achilles::CovarianceEigenShifts(kCovariance, 1e-10);
        CHECK(shifts.sigmas.size() == 5);
        CheckMatrix(FromShifts(shifts, 5), kCovariance);
    }

    SECTION("Null directions of a singular covariance are dropped") {
        // Rank two covariance in four dimensions, e.g. parameters fixed by sum rules
        const auto cov = MultiplyTransposed({{1, 0.5}, {0.2, -1}, {0.3, 0.3}, {1.2, 0.1}});
        auto shifts = achilles::CovarianceEigenShifts(cov, 1e-10);
        CHECK(shifts.sigmas.size() == 2);
        CheckMatrix(FromShifts(shifts, 4), cov);
        CHECK(achilles::CovarianceEigenShifts(cov, 1e-10, 1).sigmas.size() == 1);
    }

    SECTION("Not positive semi-definite") {
        CHECK_THROWS(achilles::CovarianceEigenShifts({{1, 2}, {2, 1}}, 1e-10));
    }

    SECTION("Hessian band of a linear observable is exact") {
        // X(p) = a . p, with sigma_X = sqrt(a^T C a)
        const std::vector<double> a{0.3, -1.2, 0.7, 2.0, -0.4};
        double variance = 0;
        for(size_t i = 0; i < 5; ++i)
            for(size_t j = 0; j < 5; ++j) variance += a[i] * kCovariance[i][j] * a[j];

        auto shifts = achilles::CovarianceEigenShifts(kCovariance, 1e-10);
        std::vector<double> members;
        for(size_t k = 0; k < shifts.sigmas.size(); ++k) {
            double projection = 0;
            for(size_t i = 0; i < 5; ++i) projection += a[i] * shifts.directions[k][i];
            members.push_back(shifts.sigmas[k] * projection);
            members.push_back(-shifts.sigmas[k] * projection);
        }
        auto band =
            achilles::Combination::Build(YAML::Load("Hessian"), "None")->Combine(0, members);
        CHECK(band.upper == Catch::Approx(std::sqrt(variance)));
        CHECK(band.lower == Catch::Approx(-std::sqrt(variance)));
    }
}

TEST_CASE("Covariance groups", "[Variations]") {
    SECTION("Form factor parameters") {
        achilles::CovarianceGroup group(YAML::Load(R"(
Name: ZExp
Type: Covariance
Target: {FormFactor: [AxialZExpansion/CC Params/1, AxialZExpansion/CC Params/2]}
Uncertainties: [0.1, 0.2]
Correlation: [[1, -0.5], [-0.5, 1]]
NSigma: 2
)"));
        CHECK(group.WeightNames() == std::vector<std::string>{"ZExp:central", "ZExp:ev1+",
                                                              "ZExp:ev1-", "ZExp:ev2+",
                                                              "ZExp:ev2-"});
        CHECK(group.GetCombination().Type() == "Hessian");
        CHECK(group.GetCombination().Scale() == Catch::Approx(0.5));
        const auto info = group.Info();
        CHECK(info.parameters == std::vector<std::string>{"AxialZExpansion/CC Params/1",
                                                          "AxialZExpansion/CC Params/2"});
        REQUIRE(info.member_values.size() == 4);
        CHECK(info.member_values[0] == Catch::Approx(-info.member_values[1]));
        REQUIRE(info.vectors.size() == 2);
        CHECK(info.vectors[0].first == "Sigmas");
        CHECK(info.vectors[1].second.size() == 4);
    }

    SECTION("Shifts are applied relative to the nominal of each model") {
        const auto nominal = BaseConfig();
        FormFactorVariation::Edit edit;
        edit.parameters = {{"AxialZExpansion/CC Params/1", 0.1, true},
                           {"AxialZExpansion/CC Params/2", -0.05, true}};
        const auto edited = FormFactorVariation::EditedConfig(nominal, edit);
        CHECK(FormFactorVariation::Lookup(edited, "AxialZExpansion/CC Params/1").as<double>() ==
              Catch::Approx(2.40));
        CHECK(FormFactorVariation::Lookup(edited, "AxialZExpansion/CC Params/2").as<double>() ==
              Catch::Approx(-0.65));
        edit.parameters = {{"AxialZExpansion/CC Params/1", 0, true}};
        CHECK(FormFactorVariation::EditedConfig(nominal, edit).IsNull());
    }

    SECTION("Flux bins") {
        achilles::CovarianceGroup group(YAML::Load(R"(
Name: Flux
Type: Covariance
Target: {Flux: {Bins: [0, 500, 1000, 2000]}}
Uncertainties: [0.1, 0.05, 0.08]
)"));
        CHECK(group.NWeights() == 7);
        CHECK(group.Info().parameters ==
              std::vector<std::string>{"Flux[0,500)", "Flux[500,1000)", "Flux[1000,2000)"});
    }

    SECTION("Shorthand") {
        CheckEquivalent(R"(
Covariance:
  ZExp:
    FormFactor: [AxialZExpansion/CC Params/1, AxialZExpansion/CC Params/2]
    Uncertainties: [0.1, 0.2]
    Correlation: [[1, -0.5], [-0.5, 1]]
  Flux:
    Flux: {Bins: [0, 500, 1000]}
    Uncertainties: [0.1, 0.05]
)",
                        R"(
- Name: ZExp
  Type: Covariance
  Target: {FormFactor: [AxialZExpansion/CC Params/1, AxialZExpansion/CC Params/2]}
  Uncertainties: [0.1, 0.2]
  Correlation: [[1, -0.5], [-0.5, 1]]
- Name: Flux
  Type: Covariance
  Target: {Flux: {Bins: [0, 500, 1000]}}
  Uncertainties: [0.1, 0.05]
)");
    }

    SECTION("Invalid groups") {
        auto build = [](const std::string &yaml) {
            return achilles::CovarianceGroup(YAML::Load(yaml));
        };
        CHECK_THROWS(build("{Name: X, Type: Covariance, Uncertainties: [0.1]}"));
        CHECK_THROWS(build("{Name: X, Type: Covariance, Target: {Fluxx: {Bins: [0, 1]}},"
                           " Uncertainties: [0.1]}"));
        CHECK_THROWS(build("{Name: X, Type: Covariance, Target: {Flux: {Bins: [1, 0]}},"
                           " Uncertainties: [0.1]}"));
        CHECK_THROWS(build("{Name: X, Type: Covariance, Target: {FormFactor: [A/b, A/c]},"
                           " Uncertainties: [0.1]}"));
    }
}
