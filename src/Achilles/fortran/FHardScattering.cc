#include <cstring>

#include "Achilles/HardScattering.hh"
#include "Achilles/HardScatteringFactory.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Event.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/Spinor.hh"

#include <fstream>

extern "C" {
    void InitializeOneBody(const char*, const char*, int, int);
    void CrossSectionOneBody(achilles::FourVector*, achilles::FourVector*,
       		                 achilles::FourVector*, achilles::FourVector*,
                             double, double, double, double, double, double,
                             unsigned long, unsigned long, double, double**, int*);
    void SetSpinors(achilles::FourVector*, achilles::FourVector*, achilles::FourVector*);
    void GetSpectral(int, double, double, double&);
    void HadronicCurrentOneBody(std::complex<double>, std::complex<double>, std::complex<double>,
                                std::complex<double>*);
    void CleanUp(double**, int*);
}

REGISTER_HARDSCATTERING(achilles::FQESpectral);
REGISTER_HARDSCATTERING(achilles::FQEGlobalFermiGas);

achilles::FQESpectral::FQESpectral(const YAML::Node &config, RunMode mode)
        : achilles::CQESpectral(mode) {

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

achilles::HCurrents achilles::FQESpectral::HadronicCurrents(Event &event,
                                                            const FFInfoMap &protonFF,
                                                            const FFInfoMap &neutronFF,
                                                            const FFInfoMap&) const {
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
    auto removal_energy = Constant::mN - pNucleonIn.E();
    auto free_energy = sqrt(pNucleonIn.P2() + Constant::mN2);
    auto ffVals = EvalFormFactor(-qVec.M2()/1.0_GeV/1.0_GeV);
    // TEST: Check the neutrino and anti-neutrion form factors
     //auto q2 = 100*100;
     //ffVals = EvalFormFactor(q2/1.0_GeV/1.0_GeV);
     //spdlog::debug("Form factor = {}, {}, {}", ffVals)
    auto omega = qVec.E();
    qVec.E() = qVec.E() + pNucleonIn.E() - free_energy;

    HCurrents results;
    double spectralProton = 0, spectralNeutron = 0;

    spdlog::debug("Energy = {}", pNucleonIn.E());
    GetSpectral(2212, removal_energy, pNucleonIn.P(), spectralProton);
    GetSpectral(2112, removal_energy, pNucleonIn.P(), spectralNeutron);
    spdlog::debug("Spectral function: S_p({}, {}) = {}, S_n({}, {}) = {}",
                  removal_energy, pNucleonIn.P(), spectralProton,                  
                  removal_energy, pNucleonIn.P(), spectralNeutron);

    // Setup spinors
    pNucleonIn.E() = free_energy;
    std::array<Spinor, 2> ubar, u;
    ubar[0] = UBarSpinor(-1, pNucleonOut);
    ubar[1] = UBarSpinor(1, pNucleonOut);
    u[0] = USpinor(-1, -pNucleonIn);
    u[1] = USpinor(1, -pNucleonIn);

    spdlog::debug("pNucleonIn = {}", pNucleonIn);
    spdlog::debug("pNucleonOut = {}", pNucleonOut);
    spdlog::debug("qVec = {}", qVec);

    // Calculate proton contributions
    for(const auto &formFactor : protonFF) {
        auto ffVal = CouplingsFF(ffVals, formFactor.second);
        spdlog::debug("f1p = {}, f2p = {}, fa = {}", ffVal[0], ffVal[1], ffVal[2]);

        auto current = CalculateHadronicCurrent(ubar, u, qVec, ffVal);

        for(auto &subcur : current) {
            for(auto &val : subcur) {
                val *= sqrt(spectralProton/6);
            }
            // Correct the Ward identity
            // subcur[3] = omega/qVec.P()*subcur[0];
        }

        results[0][formFactor.first] = current;
    }

    // Calculate neutron contributions
    for(const auto &formFactor : neutronFF) {
        auto ffVal = CouplingsFF(ffVals, formFactor.second);
        spdlog::debug("f1n = {}, f2n = {}, fa = {}", ffVal[0], ffVal[1], ffVal[2]);

        auto current = CalculateHadronicCurrent(ubar, u, qVec, ffVal);
        for(auto &subcur : current) {
            for(auto &val : subcur) {
                val *= sqrt(spectralNeutron/6);
            }
            // Correct the Ward identity
            // subcur[3] = omega/qVec.P()*subcur[0];
        }

        results[1][formFactor.first] = current;
    }

    return results;
}

void achilles::FQESpectral::CrossSection(Event &event) const {
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

    CrossSectionOneBody(&pNucleonIn, &pNucleonOut,&pLeptonIn, &pLeptonOut, pNucleonIn.E(), pNucleonIn.P(), w, qval, theta,
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

achilles::FQEGlobalFermiGas::FQEGlobalFermiGas(const YAML::Node &config, RunMode mode)
        : achilles::QEGlobalFermiGas(mode) {

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
void achilles::FQEGlobalFermiGas::CrossSection(Event &event) const {
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

    CrossSectionOneBody(&pNucleonIn, &pNucleonOut,&pLeptonIn, &pLeptonOut, pNucleonIn.E(), pNucleonIn.P(), w, qval, theta,
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

