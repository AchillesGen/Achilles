#include <cstring>

#include "nuchic/HardScattering.hh"
#include "nuchic/HardScatteringFactory.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Event.hh"
#include "nuchic/Constants.hh"

extern "C" {
    void InitializeOneBody(const char*, const char*, int, int);
    void CrossSectionOneBody(nuchic::FourVector*, nuchic::FourVector*,
                             double, double, double, double, double, double,
                             unsigned long, unsigned long, double, double**, int*);

    void InitializeResonance(const char*, const char*, int, int);
    void CrossSectionResonance(nuchic::FourVector*, nuchic::FourVector*, nuchic::FourVector*,
                               double, double, double, double, double, unsigned long, 
                               unsigned long, double, double**, double**, int*);
}

REGISTER_HARDSCATTERING(nuchic::FQESpectral);
REGISTER_HARDSCATTERING(nuchic::FQEGlobalFermiGas);
REGISTER_HARDSCATTERING(nuchic::FRSSpectral);

nuchic::FQESpectral::FQESpectral(const YAML::Node &config, RunMode mode)
        : nuchic::QESpectral(mode) {

    spdlog::trace("Initializing quasielastic spectral function model");
    auto spectralP = config["SpectralP"].as<std::string>();
    auto spectralN = config["SpectralN"].as<std::string>();
    auto iform = config["iform"].as<int>();

    size_t lenP = spectralP.size();
    size_t lenN = spectralN.size();

    auto cnameP = std::unique_ptr<char>(new char[lenP]);
    auto cnameN = std::unique_ptr<char>(new char[lenN]);
    strcpy(cnameP.get(), spectralP.c_str());
    strcpy(cnameN.get(), spectralN.c_str());

    InitializeOneBody(cnameP.get(), cnameN.get(), 0, iform);
    spdlog::trace("Finished initializing quasielastic spectral function model");
}

void nuchic::FQESpectral::CrossSection(Event &event) const {
    auto pLeptonIn = event.PhaseSpace().momentum[0];
    auto pLeptonOut = event.PhaseSpace().momentum[1];
    auto pNucleonIn = event.PhaseSpace().momentum[2];
    auto pNucleonOut = event.PhaseSpace().momentum[3];

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

    double *result = nullptr;
    int size{};

    CrossSectionOneBody(&pNucleonIn, &pNucleonOut, pNucleonIn.E(), pNucleonIn.P(), w, qval, theta,
            ee, event.CurrentNucleus() -> NProtons(), event.CurrentNucleus() -> NNucleons(),
            event.CurrentNucleus() -> FermiMomentum(0), &result, &size);

    for(size_t i = 0; i < event.MatrixElements().size(); ++i) {
        if(event.CurrentNucleus() -> Nucleons()[i].ID() == PID::proton()) {
            event.MatrixElement(i).initial_state.push_back(PID::proton());
            event.MatrixElement(i).final_state.push_back(PID::proton());
            event.MatrixElement(i).weight = result[0];
        } else {
            event.MatrixElement(i).initial_state.push_back(PID::neutron());
            event.MatrixElement(i).final_state.push_back(PID::neutron());
            event.MatrixElement(i).weight = result[1];
        }
    }
}

nuchic::FQEGlobalFermiGas::FQEGlobalFermiGas(const YAML::Node &config, RunMode mode)
        : nuchic::QEGlobalFermiGas(mode) {

    spdlog::trace("Initializing quasielastic global Fermi Gas model");
    auto spectralP = config["SpectralP"].as<std::string>();
    auto spectralN = config["SpectralN"].as<std::string>();
    auto iform = config["iform"].as<int>();

    size_t lenP = spectralP.size();
    size_t lenN = spectralN.size();

    auto cnameP = std::unique_ptr<char>(new char[lenP]);
    auto cnameN = std::unique_ptr<char>(new char[lenN]);
    strcpy(cnameP.get(), spectralP.c_str());
    strcpy(cnameN.get(), spectralN.c_str());

    InitializeOneBody(cnameP.get(), cnameN.get(), 1, iform);
    spdlog::trace("Finished initializing quasielastic global Fermi Gas model");
}

void nuchic::FQEGlobalFermiGas::CrossSection(Event &event) const {
    auto pLeptonIn = event.PhaseSpace().momentum[0];
    auto pLeptonOut = event.PhaseSpace().momentum[1];
    auto pNucleonIn = event.PhaseSpace().momentum[2];
    auto pNucleonOut = event.PhaseSpace().momentum[3];

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

    double *result = nullptr;
    int size{};

    CrossSectionOneBody(&pNucleonIn, &pNucleonOut, pNucleonIn.E(), pNucleonIn.P(), w, qval, theta,
            ee, event.CurrentNucleus() -> NProtons(), event.CurrentNucleus() -> NNucleons(),
            event.CurrentNucleus() -> FermiMomentum(0), &result, &size);

    for(size_t i = 0; i < event.MatrixElements().size(); ++i) {
        if(event.CurrentNucleus() -> Nucleons()[i].ID() == PID::proton()) {
            event.MatrixElement(i).initial_state.push_back(PID::proton());
            event.MatrixElement(i).final_state.push_back(PID::proton());
            event.MatrixElement(i).weight = result[0];
        } else {
            event.MatrixElement(i).initial_state.push_back(PID::neutron());
            event.MatrixElement(i).final_state.push_back(PID::neutron());
            event.MatrixElement(i).weight = result[1];
        }
    }
}

nuchic::FRSSpectral::FRSSpectral(const YAML::Node &config, RunMode mode)
        : nuchic::RSSpectral(mode) {

    spdlog::trace("Initializing resonance spectral function model");
    auto spectralP = config["SpectralP"].as<std::string>();
    auto spectralN = config["SpectralN"].as<std::string>();
    auto iform = config["iform"].as<int>();

    size_t lenP = spectralP.size();
    size_t lenN = spectralN.size();

    auto cnameP = std::unique_ptr<char>(new char[lenP]);
    auto cnameN = std::unique_ptr<char>(new char[lenN]);
    strcpy(cnameP.get(), spectralP.c_str());
    strcpy(cnameN.get(), spectralN.c_str());

    InitializeResonance(cnameP.get(), cnameN.get(), 0, iform);
    spdlog::trace("Finished initializing resonance spectral function model");
}

void nuchic::FRSSpectral::CrossSection(Event &event) const {
    auto pLeptonIn = event.PhaseSpace().momentum[0];
    auto pLeptonOut = event.PhaseSpace().momentum[1];
    auto pNucleonIn = event.PhaseSpace().momentum[2];
    auto pNucleonOut = event.PhaseSpace().momentum[3];
    auto pPionOut = event.PhaseSpace().momentum[4];

    auto qVec = pLeptonIn - pLeptonOut;
    auto rotMat = qVec.AlignZ();
    qVec = qVec.Rotate(rotMat);
    pLeptonIn = pLeptonIn.Rotate(rotMat);
    pLeptonOut = pLeptonOut.Rotate(rotMat);
    pNucleonIn = pNucleonIn.Rotate(rotMat);
    pNucleonOut = pNucleonOut.Rotate(rotMat);
    pPionOut = pPionOut.Rotate(rotMat);

    double ee = pLeptonIn.E();
    double theta = pLeptonIn.Angle(pLeptonOut);
    double w = qVec.E();
    double qval = qVec.P();

    double *result_p = nullptr;
    double *result_n = nullptr;
    int size{};

    CrossSectionResonance(&pNucleonIn, &pNucleonOut, &pPionOut, pNucleonIn.E(), w, qval, theta,
            ee, event.CurrentNucleus() -> NProtons(), event.CurrentNucleus() -> NNucleons(),
            event.CurrentNucleus() -> FermiMomentum(0), &result_p, &result_n, &size);

    // Channels:
    // 0: n + pi-
    // 1: n + pi0
    // 2: n + pi+
    // 3: p + pi-
    // 4: p + pi0
    // 5: p + pi+

    constexpr std::array<PID, 2> outPart1{PID::neutron(), PID::proton()};
    constexpr std::array<PID, 3> outPart2{PID::pionm(), PID::pion0(), PID::pionp()};
    for(size_t i = 0; i < event.CurrentNucleus() -> NNucleons(); ++i) {
        for(size_t j = 0; j < NChannels(); ++j) {
            if(event.CurrentNucleus() -> Nucleons()[i].ID() == PID::proton()) {
                event.MatrixElement(j+NChannels()*i).initial_state.push_back(PID::proton());
                event.MatrixElement(j+NChannels()*i).final_state.push_back(outPart1[j/3]);
                event.MatrixElement(j+NChannels()*i).final_state.push_back(outPart2[j%3]);
                event.MatrixElement(j+NChannels()*i).weight = result_p[j];
            } else {
                event.MatrixElement(j+NChannels()*i).initial_state.push_back(PID::neutron());
                event.MatrixElement(j+NChannels()*i).final_state.push_back(outPart1[j/3]);
                event.MatrixElement(j+NChannels()*i).final_state.push_back(outPart2[j%3]);
                event.MatrixElement(j+NChannels()*i).weight = result_n[j];
            }
        }
    }
}
