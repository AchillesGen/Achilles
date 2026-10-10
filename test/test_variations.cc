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
    parameter.parameter = "AxialDipole/MA";
    parameter.value = 0.9;
    auto edited = FormFactorVariation::EditedConfig(nominal, parameter);
    CHECK(edited["AxialDipole"]["MA"].as<double>() == 0.9);

    FormFactorVariation::Edit overrides;
    overrides.overrides = YAML::Load("{vector: BBBA, AxialDipole: {gan1: 1.3}}");
    edited = FormFactorVariation::EditedConfig(nominal, overrides);
    CHECK(edited["vector"].as<std::string>() == "BBBA");
    CHECK(edited["AxialDipole"]["MA"].as<double>() == 1.0);

    FormFactorVariation::Edit both = overrides;
    both.parameter = "AxialDipole/MA";
    both.value = 1.2;
    edited = FormFactorVariation::EditedConfig(nominal, both);
    CHECK(edited["AxialDipole"]["MA"].as<double>() == 1.2);
    CHECK(edited["AxialDipole"]["gan1"].as<double>() == 1.3);

    // Each edit starts from the untouched nominal
    CHECK(YAML::Dump(nominal) == reference);

    // A parameter set to its generation value is recognised as no change
    FormFactorVariation::Edit same;
    same.parameter = "AxialDipole/MA";
    same.value = 1.0;
    CHECK(FormFactorVariation::EditedConfig(nominal, same).IsNull());
}
