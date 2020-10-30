#include <cstring>
#include <iostream>


#include "nuchic/HardScattering.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"

extern "C" {
    void InitializeOneBody(const char*, const char*, int, int, int, double, int);
    double CrossSectionOneBody(int, nuchic::FourVector*, nuchic::FourVector*,
                               double, double, double, double, double, double);
}

nuchic::FQESpectral::FQESpectral(const YAML::Node &config, Beam beam,
                                 std::shared_ptr<Nucleus> nuc, HardScatteringMode mode)
    : nuchic::QESpectral(beam, nuc, mode) {

    auto spectralP = config["SpectralP"].as<std::string>();
    auto spectralN = config["SpectralN"].as<std::string>();
    auto fermiGas = config["FermiGas"].as<int>();
    auto iform = config["iform"].as<int>();
    auto nZ = nuc -> NProtons();
    auto nA = nuc -> NNeutrons();
    auto fermiMom = nuc -> FermiMomentum(0);

    size_t lenP = spectralP.size();
    size_t lenN = spectralN.size();

    auto cnameP = std::unique_ptr<char>(new char[lenP]);
    auto cnameN = std::unique_ptr<char>(new char[lenN]);
    strcpy(cnameP.get(), spectralP.c_str());
    strcpy(cnameN.get(), spectralN.c_str());

    InitializeOneBody(cnameP.get(), cnameN.get(), fermiGas, nZ, nA, fermiMom, iform);
}

double nuchic::FQESpectral::CrossSection(const Particles &particles) const {
    int inucleon = particles[2].ID() == nuchic::PID::proton() ? 1 : 2;

    auto pLeptonIn = particles[0].Momentum();
    auto pLeptonOut = particles[1].Momentum();
    auto pNucleonIn = particles[2].Momentum();
    auto pNucleonOut = particles[3].Momentum();

    auto qVec = pLeptonIn - pLeptonOut;
    auto rotMat = qVec.AlignZ();
    qVec = qVec.Rotate(rotMat);
    pLeptonIn = pLeptonIn.Rotate(rotMat);
    pLeptonOut = pLeptonOut.Rotate(rotMat);
    pNucleonIn = pNucleonIn.Rotate(rotMat);
    pNucleonOut = pNucleonOut.Rotate(rotMat);

    double ee = pLeptonIn.E();
    double theta = pLeptonIn.Angle(pLeptonOut);
    double w = qVec.E();
    double qval = qVec.P();

    return CrossSectionOneBody(inucleon, &pNucleonIn, &pNucleonOut, pNucleonIn.E(), pNucleonIn.P(), w, qval, theta, ee);
}



nuchic::FQEGlobalFermiGas::FQEGlobalFermiGas(const YAML::Node &config, Beam beam,
                                 std::shared_ptr<Nucleus> nuc, HardScatteringMode mode)
    : nuchic::QEGlobalFermiGas(beam, nuc, mode) {

    auto spectralP = config["SpectralP"].as<std::string>();
    auto spectralN = config["SpectralN"].as<std::string>();
    auto fermiGas = config["FermiGas"].as<int>();
    auto iform = config["iform"].as<int>();
    auto nZ = nuc -> NProtons();
    auto nA = nuc -> NNeutrons();
    auto fermiMom = nuc -> FermiMomentum(0);

    size_t lenP = spectralP.size();
    size_t lenN = spectralN.size();

    auto cnameP = std::unique_ptr<char>(new char[lenP]);
    auto cnameN = std::unique_ptr<char>(new char[lenN]);
    strcpy(cnameP.get(), spectralP.c_str());
    strcpy(cnameN.get(), spectralN.c_str());

    InitializeOneBody(cnameP.get(), cnameN.get(), fermiGas, nZ, nA, fermiMom, iform);
}

double nuchic::FQEGlobalFermiGas::CrossSection(const Particles &particles) const {
    int inucleon = particles[2].ID() == nuchic::PID::proton() ? 1 : 2;

    auto pLeptonIn = particles[0].Momentum();
    auto pLeptonOut = particles[1].Momentum();
    auto pNucleonIn = particles[2].Momentum();
    auto pNucleonOut = particles[3].Momentum();

    auto qVec = pLeptonIn - pLeptonOut;
    auto rotMat = qVec.AlignZ();
    qVec = qVec.Rotate(rotMat);
    pLeptonIn = pLeptonIn.Rotate(rotMat);
    pLeptonOut = pLeptonOut.Rotate(rotMat);
    pNucleonIn = pNucleonIn.Rotate(rotMat);
    pNucleonOut = pNucleonOut.Rotate(rotMat);

    double ee = pLeptonIn.E();
    double theta = pLeptonIn.Angle(pLeptonOut);
    double w = qVec.E();
    double qval = qVec.P();

    return CrossSectionOneBody(inucleon, &pNucleonIn, &pNucleonOut, pNucleonIn.E(), pNucleonIn.P(), w, qval, theta, ee);
}

