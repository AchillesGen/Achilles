#include "Achilles/FormFactor.hh"
#include "Achilles/Constants.hh"

#include "fmt/format.h"
#include "yaml-cpp/yaml.h"


achilles::FormFactor::Values achilles::FormFactor::operator()(double Q2) const {
    Values results;
    if(vector) vector -> Evaluate(Q2, results);
    if(axial) axial -> Evaluate(Q2, results);
    if(coherent) coherent -> Evaluate(Q2, results);
    return results;
}

void achilles::FormFactorImpl::Fill(double tau, FormFactor::Values &result) const {
    result.F1p = (result.Gep + tau*result.Gmp)/(1+tau);
    result.F1n = (result.Gen + tau*result.Gmn)/(1+tau);
    result.F2p = (result.Gmp - result.Gep)/(1+tau);
    result.F2n = (result.Gmn - result.Gen)/(1+tau);
}

// Vector Dipole Form Factor
achilles::VectorDipole::VectorDipole(const YAML::Node &config) {
    lambda = config["lambda"].as<double>();
    muP = config["Mu Proton"].as<double>();
    muN = config["Mu Neutron"].as<double>();
}

std::unique_ptr<achilles::FormFactorImpl> achilles::VectorDipole::Construct(achilles::FFType type,
                                                                        const YAML::Node &node) {
    Validate<VectorDipole>(type);
    return std::make_unique<VectorDipole>(node);
}

void achilles::VectorDipole::Evaluate(double Q2, FormFactor::Values &result) const {
    result.Gep = 1.0/pow(1.0+Q2/lambda/lambda, 2);
    result.Gen = -muN*Q2*result.Gep/(1+5.6*Q2/pow(Constant::mp/1_GeV, 2))/(4*pow(Constant::mp/1_GeV, 2));
    result.Gmp = muP*result.Gep;
    result.Gmn = muN*result.Gep;

    double tau = Q2/4/pow(Constant::mp/1_GeV, 2);
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
    result.FA = -gan1/pow(1.0+Q2/MA/MA, 2);
    result.FAs = gans/pow(1.0+Q2/MA/MA, 2);
}

// Vector Dummy Form Factor
achilles::VectorDummy::VectorDummy(const YAML::Node &config) {
    f1p = config["f1p"].as<double>();
    f1n = config["f1n"].as<double>();
    f2p = config["f2p"].as<double>();
    f2n = config["f2n"].as<double>();
}

std::unique_ptr<achilles::FormFactorImpl> achilles::VectorDummy::Construct(achilles::FFType type,
                                                                        const YAML::Node &node) {
    Validate<VectorDummy>(type);
    return std::make_unique<VectorDummy>(node);
}

void achilles::VectorDummy::Evaluate(double Q2, FormFactor::Values &result) const {
    result.F1p = f1p;
    result.F1n = f1n;
    result.F2p = f2p;
    result.F2n = f2n;

    double tau = Q2/4/pow(Constant::mp/1_GeV, 2);
    result.Gep = result.F1p - tau*result.F2p;
    result.Gmp = result.F1p + result.F2p;
    result.Gen = result.F1n - tau*result.F2n;
    result.Gmn = result.F1n + result.F2n;
}

// Axial-Vector Dummy Form Factor
achilles::AxialDummy::AxialDummy(const YAML::Node &config) {
    fa = config["fa"].as<double>();
    fas = config["fas"].as<double>();
}

std::unique_ptr<achilles::FormFactorImpl> achilles::AxialDummy::Construct(achilles::FFType type,
                                                                           const YAML::Node &node) {
    Validate<AxialDummy>(type);
    return std::make_unique<AxialDummy>(node);
}

void achilles::AxialDummy::Evaluate(double, FormFactor::Values &result) const {

    result.FA = fa;
    result.FAs = fas;
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
    double tau = Q2/4/pow(Constant::mp/1_GeV, 2);

    result.Gep = Parameterization(termsEp, tau);
    result.Gen = 1.0/pow(1+Q2/lambdasq, 2)*termsEn[0]*tau/(1 + termsEn[1]*tau);
    result.Gmp = muP*Parameterization(termsMp, tau);
    result.Gmn = muN*Parameterization(termsMn, tau);

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
    double tau = Q2/4/pow(Constant::mp/1_GeV, 2);

    result.Gep = Numerator(numEp, tau)/Denominator(denEp, tau);
    result.Gen = Numerator(numEn, tau)/Denominator(denEn, tau);
    result.Gmp = muP*Numerator(numMp, tau)/Denominator(denMp, tau);
    result.Gmn = muN*Numerator(numMn, tau)/Denominator(denMn, tau);

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

std::unique_ptr<achilles::FormFactorImpl> achilles::ArringtonHill::Construct(achilles::FFType type,
                                                                             const YAML::Node &node) {
    Validate<ArringtonHill>(type);
    return std::make_unique<ArringtonHill>(node);
}

void achilles::ArringtonHill::Evaluate(double Q2, FormFactor::Values &result) const {
    double tau = Q2/4/pow(Constant::mp/1_GeV, 2);

    double z = (sqrt(tcut+Q2)-sqrt(tcut-t0))/(sqrt(tcut+Q2)+sqrt(tcut-t0));

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
    const double R = 1.2*std::cbrt(A);
    r = sqrt(R*R - 5*s*s);
}

std::unique_ptr<achilles::FormFactorImpl> achilles::HelmFormFactor::Construct(achilles::FFType type,
                                                                              const YAML::Node &node) {
    Validate<HelmFormFactor>(type);
    return std::make_unique<HelmFormFactor>(node);
}

void achilles::HelmFormFactor::Evaluate(double Q2, FormFactor::Values &result) const {
    double kappa = sqrt(Q2)/Constant::HBARC;
    result.Fcoh = 3*exp(-kappa*kappa*s*s/2)*(sin(kappa*r)-kappa*r*cos(kappa*r))/pow(kappa*r, 3);
}

// Lovato Form Factor
achilles::LovatoFormFactor::LovatoFormFactor(const YAML::Node &config) {
    b = config["b"].as<double>();
    c = config["c"].as<std::array<double, 5>>();
}

std::unique_ptr<achilles::FormFactorImpl> achilles::LovatoFormFactor::Construct(achilles::FFType type,
                                                                                const YAML::Node &node) {
    Validate<LovatoFormFactor>(type);
    return std::make_unique<LovatoFormFactor>(node);
}

void achilles::LovatoFormFactor::Evaluate(double Q2, FormFactor::Values &result) const {
    spdlog::trace("LovatoFormFactor: Q2 = {}", Q2);
    double x = sqrt(Q2)/Constant::HBARC;
    result.Fcoh = exp(-0.5*pow(b*x, 2))*(c[0]+c[1]*b*x+c[2]*pow(b*x, 2)+c[3]*pow(b*x, 3)+c[4]*pow(b*x, 4))/6.0;
}
