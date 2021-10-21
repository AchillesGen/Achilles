#include <cstring>

#include "nuchic/HardScattering.hh"
#include "nuchic/HardScatteringFactory.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Event.hh"
#include "nuchic/FormFactor.hh"
#include "nuchic/Spinor.hh"

#include <fstream>

extern "C" {
    void InitializeOneBody(const char*, const char*, int, int);
    void CrossSectionOneBody(nuchic::FourVector*, nuchic::FourVector*,
                             double, double, double, double, double, double,
                             unsigned long, unsigned long, double, double**, int*);
    void SetSpinors(nuchic::FourVector*, nuchic::FourVector*, nuchic::FourVector*);
    void GetSpectral(int, double, double, double&);
    void HadronicCurrentOneBody(std::complex<double>, std::complex<double>, std::complex<double>,
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

nuchic::HCurrents nuchic::FQESpectral::HadronicCurrents(Event &event,
                                                        const FFInfoMap &protonFF,
                                                        const FFInfoMap &neutronFF) const {
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
    auto ffVals = EvalFormFactor(-qVec.M2()/1.0_GeV/1.0_GeV);

    HCurrents results;
    double spectralProton = 0, spectralNeutron = 0;

    spdlog::debug("Energy = {}", pNucleonIn.E());
    GetSpectral(2212, pNucleonIn.E(), pNucleonIn.P(), spectralProton);
    GetSpectral(2112, pNucleonIn.E(), pNucleonIn.P(), spectralNeutron);
    spdlog::trace("Spectral function: S_p({}, {}) = {}, S_n({}, {}) = {}",
                  pNucleonIn.E(), pNucleonIn.P(), spectralProton,                  
                  pNucleonIn.E(), pNucleonIn.P(), spectralNeutron);

    // Setup spinors
    pNucleonIn.E() = sqrt(pNucleonIn.P2() + Constant::mN2);
    std::array<Spinor, 2> ubar, u;
    ubar[0] = UBarSpinor(-1, pNucleonOut);
    ubar[1] = UBarSpinor(1, pNucleonOut);
    u[0] = USpinor(-1, -pNucleonIn);
    u[1] = USpinor(1, -pNucleonIn);

    spdlog::trace("ubar = {}, {}", ubar[0], ubar[1]);
    spdlog::trace("u = {}, {}", u[0], u[1]);

    // Calculate neutron contributions
    for(const auto &formFactor : protonFF) {
        std::vector<std::vector<std::complex<double>>> tmp;
        auto ffVal = CouplingsFF(ffVals, formFactor.second);
        spdlog::debug("f1p = {}, f2p = {}, fa = {}", ffVal[0], ffVal[1], ffVal[2]);

        for(size_t i = 0; i < 2; ++i) {
            for(size_t j = 0; j < 2; ++j) {
                std::vector<std::complex<double>> subcur(4);
                for(size_t mu = 0; mu < 4; ++mu) {
                    subcur[mu] = ubar[i]*(ffVal[0]*SpinMatrix::GammaMu(mu)
                                        + ffVal[2]*SpinMatrix::GammaMu(mu)*SpinMatrix::Gamma_5())*u[j];
                    double sign = 1;
                    for(size_t nu = 0; nu < 4; ++nu) {
                        subcur[mu] += ubar[i]*(ffVal[1]*SpinMatrix::SigmaMuNu(mu, nu)*sign*qVec[nu]/(2*Constant::mN))*u[j];
                        sign = -1;
                    }
                    subcur[mu] *= sqrt(spectralProton/6);
                }
                // Correct the Ward identity
                // subcur[3] = qVec.E()/qVec.P()*subcur[0];
                tmp.push_back(subcur);
            }
        }

        results[0][formFactor.first] = tmp;

    }

    // Calculate neutron contributions
    for(const auto &formFactor : neutronFF) {
        std::vector<std::vector<std::complex<double>>> tmp;
        auto ffVal = CouplingsFF(ffVals, formFactor.second);
        spdlog::debug("f1n = {}, f2n = {}, fa = {}", ffVal[0], ffVal[1], ffVal[2]);

        for(size_t i = 0; i < 2; ++i) {
            for(size_t j = 0; j < 2; ++j) {
                std::vector<std::complex<double>> subcur(4);
                for(size_t mu = 0; mu < 4; ++mu) {
                    subcur[mu] = ubar[i]*(ffVal[0]*SpinMatrix::GammaMu(mu)
                                        + ffVal[2]*SpinMatrix::GammaMu(mu)*SpinMatrix::Gamma_5())*u[j];
                    double sign = 1;
                    for(size_t nu = 0; nu < 4; ++nu) {
                        subcur[mu] += ubar[i]*(ffVal[1]*SpinMatrix::SigmaMuNu(mu, nu)*sign*qVec[nu]/(2*Constant::mN))*u[j];
                        sign = -1;
                    }

                    subcur[mu] *= sqrt(spectralNeutron/6);
                }
                // Correct the Ward identity
                // subcur[3] = qVec.E()/qVec.P()*subcur[0];
                tmp.push_back(subcur);
            }
        }

        results[1][formFactor.first] = tmp;
    }

    return results;
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

