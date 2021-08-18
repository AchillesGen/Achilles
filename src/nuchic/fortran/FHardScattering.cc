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
                                std::complex<double>*, std::complex<double>*);
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

std::pair<nuchic::Tensor, nuchic::Tensor> nuchic::FQESpectral::HadronicTensor(Event &event) const {
    auto pNucleonIn = event.Momentum().front();
    auto pNucleonOut = event.Momentum().back();
    auto qVec = event.Momentum()[1];
    for(size_t i = 2; i < event.Momentum().size()-1; ++i) {
        qVec -= event.Momentum()[i];
    }
    auto rotMat = qVec.AlignZ();
    qVec = qVec.Rotate(rotMat);
    pNucleonIn = pNucleonIn.Rotate(rotMat);
    pNucleonOut = pNucleonOut.Rotate(rotMat);
    pNucleonIn.E() = Constant::mN - pNucleonIn.E();

    Tensor result_p{}, result_n{};
    HadronicTensorOneBody(&qVec, &pNucleonIn, &pNucleonOut, result_p.data(), result_n.data());

    // Convert from fm^2 to MeV^-2
    for(size_t i = 0; i < result_p.size(); ++i) {
        result_p[i] *= pow(Constant::HBARC, 2);
        result_n[i] *= pow(Constant::HBARC, 2);
    }

    return {result_p, result_n};
}

void nuchic::FQESpectral::CrossSection(Event &event) const {
    auto pLeptonIn = event.Momentum()[1];
    auto pLeptonOut = event.Momentum()[2];
    auto pNucleonIn = event.Momentum()[0];
    auto pNucleonOut = event.Momentum()[3];

    auto qVec = pLeptonIn - pLeptonOut;
    auto rotMat = qVec.AlignZ();
    qVec = qVec.Rotate(rotMat);
    pLeptonIn = pLeptonIn.Rotate(rotMat);
    pLeptonOut = pLeptonOut.Rotate(rotMat);
    pNucleonIn = pNucleonIn.Rotate(rotMat);
    pNucleonOut = pNucleonOut.Rotate(rotMat);

    pNucleonIn.E() = Constant::mN - pNucleonIn.E();

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
            event.MatrixElement(i).inital_state[0] = PID::proton();
            event.MatrixElement(i).final_state.push_back(PID::proton());
            event.MatrixElement(i).weight = result[0];
        } else {
            event.MatrixElement(i).inital_state[0] = PID::neutron();
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

// TODO: Fix the labels to match the new ordering
void nuchic::FQEGlobalFermiGas::CrossSection(Event &event) const {
    auto pLeptonIn = event.Momentum()[0];
    auto pLeptonOut = event.Momentum()[1];
    auto pNucleonIn = event.Momentum()[2];
    auto pNucleonOut = event.Momentum()[3];

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

