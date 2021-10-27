#include "nuchic/FormFactor.hh"
#include "nuchic/Constants.hh"

#include "fmt/format.h"
#include "yaml-cpp/yaml.h"

std::unique_ptr<nuchic::FormFactor> nuchic::FormFactor::Build(const YAML::Node &node) {
    std::string formFactor = node["FormFactor"].as<std::string>();

    std::unique_ptr<FormFactor> result;
    if(formFactor == "Dipole")
        result = std::make_unique<Dipole>(Dipole(node["Dipole"]));
    else if(formFactor == "Kelly")
        result = std::make_unique<Kelly>(Kelly(node["Kelly"]));
    else if(formFactor == "BBBA")
        result = std::make_unique<BBBA>(BBBA(node["BBBA"]));
    else if(formFactor == "ArringtonHill")
        result = std::make_unique<ArringtonHill>(ArringtonHill(node["ArringtonHill"]));
    else
        throw std::runtime_error(fmt::format("Invalid Form Factor: {}", formFactor));

    return result;
}

void nuchic::FormFactor::Fill(double tau, Values &result) const {
    result.F1p = (result.Gep + tau*result.Gmp)/(1+tau);
    result.F1n = (result.Gen + tau*result.Gmn)/(1+tau);
    result.F2p = (result.Gmp - result.Gep)/(1+tau);
    result.F2n = (result.Gmn - result.Gen)/(1+tau);
}

nuchic::Dipole::Dipole(const YAML::Node &config) {
    lambda = config["lambda"].as<double>();
    MA = config["MA"].as<double>();
    muP = config["Mu Proton"].as<double>();
    muN = config["Mu Neutron"].as<double>();
    gan1 = config["gan1"].as<double>();
    gans = config["gans"].as<double>();
}

// Q2 in GeV^2
nuchic::FormFactor::Values nuchic::Dipole::operator()(double Q2) const {
    Values result{};

    result.Gep = 1.0/pow(1.0+Q2/lambda/lambda, 2);
    result.Gen = -muN*Q2*result.Gep/(1+5.6*Q2/pow(Constant::mp/1_GeV, 2))/(4*pow(Constant::mp/1_GeV, 2));
    result.Gmp = muP*result.Gep;
    result.Gmp = muN*result.Gep;
    result.FA = -gan1/pow(1.0+Q2/MA/MA, 2);
    result.FAs = -gans/pow(1.0+Q2/MA/MA, 2);

    double tau = Q2/4/pow(Constant::mp/1_GeV, 2);
    Fill(tau, result);
    return result;
}

nuchic::Kelly::Kelly(const YAML::Node &config) {
    lambdasq = config["lambdasq"].as<double>();
    MA = config["MA"].as<double>();
    muP = config["Mu Proton"].as<double>();
    muN = config["Mu Neutron"].as<double>();
    termsEp = config["Gep Params"].as<std::array<double, 4>>();
    termsEn = config["Gen Params"].as<std::array<double, 2>>();
    termsMp = config["Gmp Params"].as<std::array<double, 4>>();
    termsMn = config["Gmn Params"].as<std::array<double, 4>>();
    gan1 = config["gan1"].as<double>();
    gans = config["gans"].as<double>();
}

nuchic::FormFactor::Values nuchic::Kelly::operator()(double Q2) const {
    Values result{};

    double tau = Q2/4/pow(Constant::mp/1_GeV, 2);

    result.Gep = Parameterization(termsEp, tau);
    result.Gen = 1.0/pow(1+Q2/lambdasq, 2)*termsEn[0]*tau/(1 + termsEn[1]*tau);
    result.Gmp = muP*Parameterization(termsMp, tau);
    result.Gmn = muN*Parameterization(termsMn, tau);
    result.FA = -gan1/pow(1.0+Q2/MA/MA, 2);
    result.FAs = -gans/pow(1.0+Q2/MA/MA, 2);

    Fill(tau, result);
    return result;
}

nuchic::BBBA::BBBA(const YAML::Node &config) {
    muP = config["Mu Proton"].as<double>();
    muN = config["Mu Neutron"].as<double>();
    MA = config["MA"].as<double>();
    gan1 = config["gan1"].as<double>();
    gans = config["gans"].as<double>();

    numEp = config["NumeratorEp Params"].as<std::array<double, 4>>();
    denEp = config["DenominatorEp Params"].as<std::array<double, 4>>();
    numEn = config["NumeratorEn Params"].as<std::array<double, 4>>();
    denEn = config["DenominatorEn Params"].as<std::array<double, 4>>();
    numMp = config["NumeratorMp Params"].as<std::array<double, 4>>();
    denMp = config["DenominatorMp Params"].as<std::array<double, 4>>();
    numMn = config["NumeratorMn Params"].as<std::array<double, 4>>();
    denMn = config["DenominatorMn Params"].as<std::array<double, 4>>();
}

nuchic::FormFactor::Values nuchic::BBBA::operator()(double Q2) const {
    Values result{};

    double tau = Q2/4/pow(Constant::mp/1_GeV, 2);

    result.Gep = Numerator(numEp, tau)/Denominator(denEp, tau);
    result.Gen = Numerator(numEn, tau)/Denominator(denEn, tau);
    result.Gmp = muP*Numerator(numMp, tau)/Denominator(denMp, tau);
    result.Gmn = muN*Numerator(numMn, tau)/Denominator(denMn, tau);
    result.FA = -gan1/pow(1.0+Q2/MA/MA, 2);
    result.FAs = -gans/pow(1.0+Q2/MA/MA, 2);

    Fill(tau, result);
    return result;
}

nuchic::ArringtonHill::ArringtonHill(const YAML::Node &config) {
    muP = config["Mu Proton"].as<double>();
    muN = config["Mu Neutron"].as<double>();
    MA = config["MA"].as<double>();
    gan1 = config["gan1"].as<double>();
    gans = config["gans"].as<double>();

    tcut = config["tcut"].as<double>();
    t0 = config["t0"].as<double>();

    epParams = config["Gep Params"].as<std::array<double, 13>>();
    enParams = config["Gen Params"].as<std::array<double, 13>>();
    mpParams = config["Gmp Params"].as<std::array<double, 13>>();
    mnParams = config["Gmn Params"].as<std::array<double, 13>>();
}

nuchic::FormFactor::Values nuchic::ArringtonHill::operator()(double Q2) const {
    Values result{};

    double tau = Q2/4/pow(Constant::mp/1_GeV, 2);

    double z = (sqrt(tcut+Q2)-sqrt(tcut-t0))/(sqrt(tcut+Q2)+sqrt(tcut-t0));

    result.Gep = ZExpand(epParams, z);
    result.Gen = ZExpand(enParams, z);
    result.Gmp = ZExpand(mpParams, z);
    result.Gmn = ZExpand(mnParams, z);
    result.FA = -gan1/pow(1.0+Q2/MA/MA, 2);
    result.FAs = -gans/pow(1.0+Q2/MA/MA, 2);

    Fill(tau, result);
    return result;
}
