#include "Achilles/FormFactor.hh"
#include "Achilles/Constants.hh"

#include "Achilles/ParticleInfo.hh"
#include "fmt/format.h"
#include "yaml-cpp/yaml.h"

achilles::FormFactor::Values achilles::FormFactor::operator()(double Q2) const {
    Values results;
    if(vector) vector->Evaluate(Q2, results);
    if(axial) axial->Evaluate(Q2, results);
    if(coherent) coherent->Evaluate(Q2, results);
    if(resonancevector) resonancevector->Evaluate(Q2, results);
    if(resonanceaxial) resonanceaxial->Evaluate(Q2, results);
    if(mecvector) mecvector->Evaluate(Q2, results);
    if(mecaxial) mecaxial->Evaluate(Q2, results);
    if(hyperon) hyperon->Evaluate(Q2, results);
    return results;
}

void achilles::FormFactorImpl::Fill(double tau, FormFactor::Values &result) const {
    result.F1p = (result.Gep + tau * result.Gmp) / (1 + tau);
    result.F1n = (result.Gen + tau * result.Gmn) / (1 + tau);
    result.F2p = (result.Gmp - result.Gep) / (1 + tau);
    result.F2n = (result.Gmn - result.Gen) / (1 + tau);
}

// Vector Dipole Form Factor
achilles::VectorDipole::VectorDipole(const YAML::Node &config) {
    lambda = config["lambda"].as<double>();
    muP = config["Mu Proton"].as<double>();
    muN = config["Mu Neutron"].as<double>();
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::VectorDipole::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<VectorDipole>(type);
    return std::make_unique<VectorDipole>(node);
}

void achilles::VectorDipole::Evaluate(double Q2, FormFactor::Values &result) const {
    result.Gep = 1.0 / pow(1.0 + Q2 / lambda / lambda, 2);
    result.Gen = -muN * Q2 * result.Gep / (1 + 5.6 * Q2 / pow(Constant::mp / 1_GeV, 2)) /
                 (4 * pow(Constant::mp / 1_GeV, 2));
    result.Gmp = muP * result.Gep;
    result.Gmn = muN * result.Gep;

    double tau = Q2 / 4 / pow(Constant::mp / 1_GeV, 2);
    Fill(tau, result);
}

// Axial-Vector Dipole Form Factor
achilles::AxialDipole::AxialDipole(const YAML::Node &config) {
    MA = config["MA"].as<double>();
    gan1 = config["gan1"].as<double>();
    gans = config["gans"].as<double>();
}

std::unique_ptr<achilles::FormFactorImpl> achilles::AxialDipole::Construct(achilles::FFType type,
                                                                           const YAML::Node &node) {
    Validate<AxialDipole>(type);
    return std::make_unique<AxialDipole>(node);
}

void achilles::AxialDipole::Evaluate(double Q2, FormFactor::Values &result) const {
    result.FA = -gan1 / pow(1.0 + Q2 / MA / MA, 2);
    result.FAs = -gans / pow(1.0 + Q2 / MA / MA, 2);
    result.FAP =
        2.0 * Constant::mN2 / (Q2 * 1_GeV * 1_GeV + pow(ParticleInfo(211).Mass(), 2)) * result.FA;
}

// Axial Z-Expansion Form Factor
achilles::AxialZExpansion::AxialZExpansion(const YAML::Node &config) {
    tcut = config["tcut"].as<double>();
    t0 = config["t0"].as<double>();
    cc_params = config["CC Params"].as<std::vector<double>>();
    strange_params = config["Strange Params"].as<std::vector<double>>();
    spdlog::debug("Parameters for AxialZExpansion");
    for(size_t i = 0; i < cc_params.size(); i++) {
        spdlog::debug("cc_params[{0}] = {1}", i, cc_params[i]);
    }
    for(size_t i = 0; i < strange_params.size(); i++) {
        spdlog::debug("strange_params[{0}] = {1}", i, strange_params[i]);
    }
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::AxialZExpansion::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<AxialZExpansion>(type);
    return std::make_unique<AxialZExpansion>(node);
}

void achilles::AxialZExpansion::Evaluate(double Q2, FormFactor::Values &result) const {
    double z = (sqrt(tcut + Q2) - sqrt(tcut - t0)) / (sqrt(tcut + Q2) + sqrt(tcut - t0));
    // Charged-current axial form factor
    result.FA = ZExpand(cc_params, z);
    // Strange-quark axial form factor
    result.FAs = ZExpand(strange_params, z);
    // Pseudo-axial form factor
    result.FAP =
        2.0 * Constant::mN2 / (Q2 * 1_GeV * 1_GeV + pow(ParticleInfo(211).Mass(), 2)) * result.FA;
}

// Kelly Form Factor
achilles::Kelly::Kelly(const YAML::Node &config) {
    lambdasq = config["lambdasq"].as<double>();
    muP = config["Mu Proton"].as<double>();
    muN = config["Mu Neutron"].as<double>();
    termsEp = config["Gep Params"].as<std::array<double, 4>>();
    termsEn = config["Gen Params"].as<std::array<double, 2>>();
    termsMp = config["Gmp Params"].as<std::array<double, 4>>();
    termsMn = config["Gmn Params"].as<std::array<double, 4>>();
}

std::unique_ptr<achilles::FormFactorImpl> achilles::Kelly::Construct(achilles::FFType type,
                                                                     const YAML::Node &node) {
    Validate<Kelly>(type);
    return std::make_unique<Kelly>(node);
}

void achilles::Kelly::Evaluate(double Q2, FormFactor::Values &result) const {
    double tau = Q2 / 4 / pow(Constant::mp / 1_GeV, 2);

    result.Gep = Parameterization(termsEp, tau);
    result.Gen = 1.0 / pow(1 + Q2 / lambdasq, 2) * termsEn[0] * tau / (1 + termsEn[1] * tau);
    result.Gmp = muP * Parameterization(termsMp, tau);
    result.Gmn = muN * Parameterization(termsMn, tau);

    Fill(tau, result);
}

// BBBA Form Factor
achilles::BBBA::BBBA(const YAML::Node &config) {
    muP = config["Mu Proton"].as<double>();
    muN = config["Mu Neutron"].as<double>();

    numEp = config["NumeratorEp Params"].as<std::array<double, 4>>();
    denEp = config["DenominatorEp Params"].as<std::array<double, 4>>();
    numEn = config["NumeratorEn Params"].as<std::array<double, 4>>();
    denEn = config["DenominatorEn Params"].as<std::array<double, 4>>();
    numMp = config["NumeratorMp Params"].as<std::array<double, 4>>();
    denMp = config["DenominatorMp Params"].as<std::array<double, 4>>();
    numMn = config["NumeratorMn Params"].as<std::array<double, 4>>();
    denMn = config["DenominatorMn Params"].as<std::array<double, 4>>();
}

std::unique_ptr<achilles::FormFactorImpl> achilles::BBBA::Construct(achilles::FFType type,
                                                                    const YAML::Node &node) {
    Validate<BBBA>(type);
    return std::make_unique<BBBA>(node);
}

void achilles::BBBA::Evaluate(double Q2, FormFactor::Values &result) const {
    double tau = Q2 / 4 / pow(Constant::mp / 1_GeV, 2);

    result.Gep = Numerator(numEp, tau) / Denominator(denEp, tau);
    result.Gen = Numerator(numEn, tau) / Denominator(denEn, tau);
    result.Gmp = muP * Numerator(numMp, tau) / Denominator(denMp, tau);
    result.Gmn = muN * Numerator(numMn, tau) / Denominator(denMn, tau);

    Fill(tau, result);
}

// Arrington-Hill Form Factor
achilles::ArringtonHill::ArringtonHill(const YAML::Node &config) {
    muP = config["Mu Proton"].as<double>();
    muN = config["Mu Neutron"].as<double>();

    tcut = config["tcut"].as<double>();
    t0 = config["t0"].as<double>();

    epParams = config["Gep Params"].as<std::array<double, 13>>();
    enParams = config["Gen Params"].as<std::array<double, 13>>();
    mpParams = config["Gmp Params"].as<std::array<double, 13>>();
    mnParams = config["Gmn Params"].as<std::array<double, 13>>();
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::ArringtonHill::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<ArringtonHill>(type);
    return std::make_unique<ArringtonHill>(node);
}

void achilles::ArringtonHill::Evaluate(double Q2, FormFactor::Values &result) const {
    double tau = Q2 / 4 / pow(Constant::mp / 1_GeV, 2);

    double z = (sqrt(tcut + Q2) - sqrt(tcut - t0)) / (sqrt(tcut + Q2) + sqrt(tcut - t0));

    result.Gep = ZExpand(epParams, z);
    result.Gen = ZExpand(enParams, z);
    result.Gmp = ZExpand(mpParams, z);
    result.Gmn = ZExpand(mnParams, z);

    Fill(tau, result);
}

// Helm Form Factor
achilles::HelmFormFactor::HelmFormFactor(const YAML::Node &config) {
    s = config["s"].as<double>();
    auto A = config["A"].as<double>();
    const double R = 1.2 * std::cbrt(A);
    r = sqrt(R * R - 5 * s * s);
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::HelmFormFactor::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<HelmFormFactor>(type);
    return std::make_unique<HelmFormFactor>(node);
}

void achilles::HelmFormFactor::Evaluate(double Q2, FormFactor::Values &result) const {
    double kappa = sqrt(Q2) / Constant::HBARC;
    result.Fcoh = 3 * exp(-kappa * kappa * s * s / 2) *
                  (sin(kappa * r) - kappa * r * cos(kappa * r)) / pow(kappa * r, 3);
}

// Klein-Nystrand (KN) Form Factor 2007.03658 (Eq. 26)
achilles::KNFormFactor::KNFormFactor(const YAML::Node &config) {
    r0 = config["r0"].as<double>(); // default = 3.427 fm
    ak = config["ak"].as<double>(); // default = 0.7 fm
    if(config["Adapted"] && config["Adapted"].as<bool>()) {
        RA = sqrt(5 * r0 * r0 / 3 - 10 * ak * ak);
    } else {
        auto A = config["A"].as<double>();
        RA = 1.2 * std::cbrt(A); // default = 1.2 * A^(1/3) fm
    }
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::KNFormFactor::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<KNFormFactor>(type);
    return std::make_unique<KNFormFactor>(node);
}

void achilles::KNFormFactor::Evaluate(double Q2, FormFactor::Values &result) const {
    double kappa = sqrt(Q2) / Constant::HBARC;
    result.Fcoh = 3 / (1 + kappa * kappa * ak * ak) *
                  (sin(kappa * RA) - kappa * RA * cos(kappa * RA)) / pow(kappa * RA, 3);
}

// Lovato Form Factor
achilles::LovatoFormFactor::LovatoFormFactor(const YAML::Node &config) {
    b = config["b"].as<double>();
    c = config["c"].as<std::array<double, 5>>();
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::LovatoFormFactor::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<LovatoFormFactor>(type);
    return std::make_unique<LovatoFormFactor>(node);
}

void achilles::LovatoFormFactor::Evaluate(double Q2, FormFactor::Values &result) const {
    spdlog::trace("LovatoFormFactor: Q2 = {}", Q2);
    double x = sqrt(Q2) / Constant::HBARC;
    result.Fcoh =
        exp(-0.5 * pow(b * x, 2)) *
        (c[0] + c[1] * b * x + c[2] * pow(b * x, 2) + c[3] * pow(b * x, 3) + c[4] * pow(b * x, 4)) /
        6.0;
}

// Resonance Dummy Form Factor
achilles::ResonanceDummyVectorFormFactor::ResonanceDummyVectorFormFactor(const YAML::Node &config) {
    resV = config["resV"].as<double>();
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::ResonanceDummyVectorFormFactor::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<ResonanceDummyVectorFormFactor>(type);
    return std::make_unique<ResonanceDummyVectorFormFactor>(node);
}

void achilles::ResonanceDummyVectorFormFactor::Evaluate(double Q2,
                                                        FormFactor::Values &result) const {
    spdlog::trace("ResonanceDummyVectorFormFactor: Q2 = {}", Q2);
    result.FresV = resV;
}

// Resonance Dummy Form Factor
achilles::ResonanceDummyAxialFormFactor::ResonanceDummyAxialFormFactor(const YAML::Node &config) {
    resA = config["resA"].as<double>();
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::ResonanceDummyAxialFormFactor::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<ResonanceDummyAxialFormFactor>(type);
    return std::make_unique<ResonanceDummyAxialFormFactor>(node);
}

void achilles::ResonanceDummyAxialFormFactor::Evaluate(double Q2,
                                                       FormFactor::Values &result) const {
    spdlog::trace("ResonanceDummyAxialFormFactor: Q2 = {}", Q2);
    result.FresA = resA;
}

// MEC Vector form factors
achilles::MECVectorFormFactor::MECVectorFormFactor(const YAML::Node &config) {
    MvSq = config["MvSq"].as<double>();
    cv3norm = config["cv3norm"].as<double>();
    cv4norm = config["cv4norm"].as<double>();
    cv5norm = config["cv5norm"].as<double>();
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::MECVectorFormFactor::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<MECVectorFormFactor>(type);
    return std::make_unique<MECVectorFormFactor>(node);
}

void achilles::MECVectorFormFactor::Evaluate(double Q2, FormFactor::Values &result) const {
    spdlog::trace("MECVectorFormFactor: Q2 = {}", Q2);
    result.FmecV3 = cv3norm / pow(1. + Q2 / MvSq, 2) / (1. + Q2 / 4. / MvSq) * sqrt(3. / 2.);
    result.FmecV4 = cv4norm / pow(1. + Q2 / MvSq, 2) / (1. + Q2 / 4. / MvSq) * sqrt(3. / 2.);
    result.FmecV5 = cv5norm / pow(1. + Q2 / MvSq, 2) / (1. + Q2 / 0.776 / MvSq) * sqrt(3. / 2.);
    result.Fpiem = result.Gep - result.Gen;
}

// MEC Axial form factors
achilles::MECAxialFormFactor::MECAxialFormFactor(const YAML::Node &config) {
    MaDeltaSq = config["MaDeltaSq"].as<double>();
    ca5norm = config["ca5norm"].as<double>();
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::MECAxialFormFactor::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<MECAxialFormFactor>(type);
    return std::make_unique<MECAxialFormFactor>(node);
}

void achilles::MECAxialFormFactor::Evaluate(double Q2, FormFactor::Values &result) const {
    spdlog::trace("MECAxialFormFactor: Q2 = {}", Q2);
    result.FmecA5 =
        ca5norm / pow(1. + Q2 / MaDeltaSq, 2) / (1. + Q2 / 3. / MaDeltaSq) * sqrt(3. / 2.);
}

// Hyperon form factors
achilles::HyperonFormFactor::HyperonFormFactor(const YAML::Node &config) {
    dummy = config["dummy"].as<double>();
}

std::unique_ptr<achilles::FormFactorImpl>
achilles::HyperonFormFactor::Construct(achilles::FFType type, const YAML::Node &node) {
    Validate<HyperonFormFactor>(type);
    return std::make_unique<HyperonFormFactor>(node);
}

void achilles::HyperonFormFactor::Evaluate(double Q2, FormFactor::Values &result) const {
    spdlog::trace("HyperonFormFactor: Q2 = {}", Q2);

    double x = 0.73;
    double lambda_ratio = Constant::mlambda / (Constant::mlambda + Constant::mn);
    double sigmam_ratio = Constant::msigmam / (Constant::msigmam + Constant::mn);
    double sigma0_ratio = Constant::msigma0 / (Constant::msigma0 + Constant::mn);

    result.F1lam = -sqrt(3. / 2.) * result.F1p;
    result.F2lam = -sqrt(3. / 2.) * lambda_ratio * result.F2p;
    result.FAlam = -sqrt(3. / 2.) * ((1. + 2. * x) / 3.) * result.FA;
    result.F1sigm = -(result.F1p + 2 * result.F1n);
    result.F2sigm = -sigmam_ratio * (result.F2p + 2. * result.F2n);
    result.FAsigm = (1. - 2. * x) * result.FA;
    result.F1sig0 = -(1. / sqrt(2.)) * (result.F1p + 2. * result.F1n);
    result.F2sig0 = -sigma0_ratio * (1. / sqrt(2.)) * (result.F2p + 2. * result.F2n);
    result.FAsig0 = (1. - 2. * x) * result.FA / sqrt(2.);
}
