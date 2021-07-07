#include <cstring>

#include "nuchic/HardScattering.hh"
#include "nuchic/HardScatteringFactory.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Event.hh"

extern "C" {
    void InitializeOneBody(const char*, const char*, int, int);
    void CrossSectionOneBody(nuchic::FourVector*, nuchic::FourVector*,
                             double, double, double, double, double, double,
                             unsigned long, unsigned long, double, double**, int*);
    void HadronicTensorOneBody(nuchic::FourVector*, nuchic::FourVector*, nuchic::FourVector*,
                                std::complex<double>*);
    void CleanUp(double**, int*);
}

REGISTER_HARDSCATTERING(nuchic::FQESpectral);
REGISTER_HARDSCATTERING(nuchic::FQEGlobalFermiGas);

nuchic::FQESpectral::FQESpectral(const YAML::Node &config, RunMode mode)
        : nuchic::QESpectral(mode) {

    spdlog::trace("Initializing quasielastic spectral function model");
    auto spectralP = config["SpectralP"].as<std::string>();
    auto spectralN = config["SpectralN"].as<std::string>();
    auto iform = config["iform"].as<int>();

    size_t lenP = spectralP.size();
    size_t lenN = spectralN.size();

    auto cnameP = std::make_unique<char[]>(lenP+1);
    auto cnameN = std::make_unique<char[]>(lenN+1);
    strcpy(cnameP.get(), spectralP.c_str());
    strcpy(cnameN.get(), spectralN.c_str());

    InitializeOneBody(cnameP.get(), cnameN.get(), 0, iform);
    spdlog::trace("Finished initializing quasielastic spectral function model");
}

nuchic::Tensor nuchic::FQESpectral::HadronicTensor(Event &event) const {
    auto pNucleonIn = event.PhaseSpace().momentum[NLeptons()];
    auto pNucleonOut = event.PhaseSpace().momentum[NLeptons()+1];
    auto qVec = event.PhaseSpace().momentum[0];
    for(size_t i = 1; i < NLeptons(); ++i) {
        qVec -= event.PhaseSpace().momentum[i];
    }
    auto rotMat = qVec.AlignZ();
    qVec = qVec.Rotate(rotMat);
    pNucleonIn = pNucleonIn.Rotate(rotMat);
    pNucleonOut = pNucleonOut.Rotate(rotMat);

    Tensor result{};
    HadronicTensorOneBody(&qVec, &pNucleonIn, &pNucleonOut, result.data());

    return result;
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
            event.MatrixElement(i).inital_state.push_back(PID::proton());
            event.MatrixElement(i).final_state.push_back(PID::proton());
            event.MatrixElement(i).weight = result[0];
        } else {
            event.MatrixElement(i).inital_state.push_back(PID::neutron());
            event.MatrixElement(i).final_state.push_back(PID::neutron());
            event.MatrixElement(i).weight = result[1];
        }
    }

    CleanUp(&result, &size);
    result = nullptr;
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
            event.MatrixElement(i).inital_state.push_back(PID::proton());
            event.MatrixElement(i).final_state.push_back(PID::proton());
            event.MatrixElement(i).weight = result[0];
        } else {
            event.MatrixElement(i).inital_state.push_back(PID::neutron());
            event.MatrixElement(i).final_state.push_back(PID::neutron());
            event.MatrixElement(i).weight = result[1];
        }
    }

    CleanUp(&result, &size);
    result = nullptr;
}

