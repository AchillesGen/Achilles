#include <iostream>
#include <fstream>
#include <map>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

#include "nuchic/Constants.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/InteractionComponent.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"
#include "nuchic/Utilities.hh"

using namespace H5;
using namespace nuchic;

REGISTER_INTERACTION_COMPONENT(GeantInteractions);
REGISTER_INTERACTION_COMPONENT(NasaInteractions);
REGISTER_INTERACTION_COMPONENT(ConstantInteractions);
REGISTER_INTERACTION_COMPONENT(PionNucleon);

const std::map<std::string, double> HZETRN = {
    {"a", 5.0_MeV},
    {"b", 0.199/sqrt(1_MeV)},
    {"c", 0.451 * pow(1_MeV, -0.258)},
    {"d", 25.0_MeV},
    {"e", 134.0_MeV},
    {"f", 1.187 * pow(1_MeV, -0.35)},
    {"g", 0.1_MeV},
    {"h", 0.282_MeV}
};

const std::map<std::string, double> PDG = {
    {"Zpp", 33.45_mb},  
    {"Zpn", 35.80_mb},  
    {"Y1pp", 42.53_mb}, 
    {"Y1pn", 40.15_mb}, 
    {"Y2pp", 33.34_mb}, 
    {"Y2pn", 30.00_mb}, 
    {"B", 0.308_mb},    
    {"s1", 1.0 * pow(1_GeV, 2)},
    {"s0", pow(5.38_GeV, 2)},
    {"n1", 0.458},
    {"n2", 0.545}
};

const std::map<std::string, double> JWN = {
    {"gamma", 52.5 * pow(1_GeV, 0.16)}, //mb
    {"alpha", 0.00369 / 1_MeV},
    {"beta", 0.00895741 * pow(1_MeV, -0.8)}
};

double nuchic::CrossSection(bool, const double&) {
    // TODO: Implement GEANT4 cross-section or something else
    throw std::domain_error("Invalid energy!");
}

double nuchic::CrossSectionLab(bool samePID, const double& pLab) noexcept {
    const double tLab = sqrt(pow(pLab, 2) + pow(Constant::mN, 2)) - Constant::mN;
    if(samePID) {
        if(pLab < 1.8_GeV) {
            if(tLab >= 25_MeV)
                return (1.0+HZETRN.at("a")/tLab) * (40+109.0*std::cos(HZETRN.at("b")*sqrt(tLab))
                        * exp(-HZETRN.at("c")*pow(tLab-HZETRN.at("d"), 0.258)));
            else
                return exp(6.51*exp(-pow(tLab/HZETRN.at("e"), 0.7)));
        } else if(pLab <= 4.7_GeV) {
            return JWN.at("gamma")/pow(pLab, 0.16);
        } else {
            double ecm2 = 2*Constant::mN*(Constant::mN+sqrt(pow(pLab, 2) + pow(Constant::mN, 2)));
            return PDG.at("Zpp") + PDG.at("B")*pow(log(ecm2/PDG.at("s0")), 2)
                + PDG.at("Y1pp")*pow(PDG.at("s1")/ecm2, PDG.at("n1"))
                - PDG.at("Y2pp")*pow(PDG.at("s1")/ecm2, PDG.at("n2"));
        }
    } else {
        if(pLab < 0.5_GeV) {
            if(tLab >= 0.1_MeV)
                return 38.0 + 12500.0*exp(-HZETRN.at("f")*pow(tLab-HZETRN.at("g"), 0.35));
            else
                return 26000 * exp(-pow(tLab/HZETRN.at("h"), 0.3));
        } else if(pLab <= 2.0_GeV) {
            return 40 + 10*cos(JWN.at("alpha")*pLab - 0.943)
                * exp(-JWN.at("beta")*pow(pLab, 0.8)+2);
        } else {
            double ecm2 = 2*Constant::mN*(Constant::mN+sqrt(pow(pLab, 2) + pow(Constant::mN, 2)));
            return PDG.at("Zpn") + PDG.at("B")*pow(log(ecm2/PDG.at("s0")), 2)
                + PDG.at("Y1pn")*pow(PDG.at("s1")/ecm2, PDG.at("n1"))
                - PDG.at("Y2pn")*pow(PDG.at("s1")/ecm2, PDG.at("n2"));
        }
    }
}

double nuchic::CrossSectionAngle(bool, const double&, const double& ran) {
    // For testing right now
    // TODO: Implement GEANT4 angular distribution or something else
    return std::acos(2*ran - 1);
}

ThreeVector nuchic::MakeMomentumAngular(bool samePID, const double& p1CM, const double& pcm,
        const std::array<double, 2>& rans) {
    double pR = p1CM;
    double pTheta = nuchic::CrossSectionAngle(samePID, pcm, rans[0]);
    double pPhi = 2*M_PI*rans[1];

    return ThreeVector(ToCartesian({pR, pTheta, pPhi}));
}

double InteractionComponent::CrossSectionLab(bool samePID, const double& pLab) const noexcept {
    const double tLab = sqrt(pow(pLab, 2) + pow(Constant::mN, 2)) - Constant::mN;
    if(samePID) {
        if(pLab < 1.8_GeV) {
            if(tLab >= 25_MeV)
                return (1.0+HZETRN.at("a")/tLab) * (40+109.0*std::cos(HZETRN.at("b")*sqrt(tLab))
                        * exp(-HZETRN.at("c")*pow(tLab-HZETRN.at("d"), 0.258)));
            else
                return exp(6.51*exp(-pow(tLab/HZETRN.at("e"), 0.7)));
        } else if(pLab <= 4.7_GeV) {
            return JWN.at("gamma")/pow(pLab, 0.16);
        } else {
            double ecm2 = 2*Constant::mN*(Constant::mN+sqrt(pow(pLab, 2) + pow(Constant::mN, 2)));
            return PDG.at("Zpp") + PDG.at("B")*pow(log(ecm2/PDG.at("s0")), 2)
                + PDG.at("Y1pp")*pow(PDG.at("s1")/ecm2, PDG.at("n1"))
                - PDG.at("Y2pp")*pow(PDG.at("s1")/ecm2, PDG.at("n2"));
        }
    } else {
        if(pLab < 0.5_GeV) {
            if(tLab >= 0.1_MeV)
                return 38.0 + 12500.0*exp(-HZETRN.at("f")*pow(tLab-HZETRN.at("g"), 0.35));
            else
                return 26000 * exp(-pow(tLab/HZETRN.at("h"), 0.3));
        } else if(pLab <= 2.0_GeV) {
            return 40 + 10*cos(JWN.at("alpha")*pLab - 0.943)
                * exp(-JWN.at("beta")*pow(pLab, 0.8)+2);
        } else {
            double ecm2 = 2*Constant::mN*(Constant::mN+sqrt(pow(pLab, 2) + pow(Constant::mN, 2)));
            return PDG.at("Zpn") + PDG.at("B")*pow(log(ecm2/PDG.at("s0")), 2)
                + PDG.at("Y1pn")*pow(PDG.at("s1")/ecm2, PDG.at("n1"))
                - PDG.at("Y2pn")*pow(PDG.at("s1")/ecm2, PDG.at("n2"));
        }
    }
}

GeantInteractions::GeantInteractions(const YAML::Node& node) {
    auto filename = node["GeantData"].as<std::string>();

    // Initialize theta vector
    constexpr double thetaMin = 0.5;
    constexpr double thetaMax = 179.5;
    constexpr size_t nTheta = 180;
    m_theta = Linspace(thetaMin, thetaMax, nTheta);

    constexpr double cdfMin = -3;
    constexpr double cdfMax = 0;
    constexpr size_t nCDF = 200;
    m_cdf = Logspace(cdfMin, cdfMax, nCDF);

    // Read in the Geant4 hdf5 file and get the np and pp groups
    spdlog::info("GeantInteractions: Loading Geant4 data from {0}.", filename);
    H5File file(filename, H5F_ACC_RDONLY);
    Group dataNP(file.openGroup("np")); 
    Group dataPP(file.openGroup("pp"));

    // Get the datasets for np and load into local variables
    LoadData(false, dataNP);

    // Get the datasets for pp and load into local variables
    LoadData(true, dataPP);
}

void GeantInteractions::LoadData(bool samePID, const Group& group) {
    // Load datasets
    DataSet pcm(group.openDataSet("pcm"));
    DataSet sigTot(group.openDataSet("sigtot"));
    DataSet sig(group.openDataSet("sig"));
  
    // Get data for center of momentum
    DataSpace pcmSpace = pcm.getSpace();
    std::array<hsize_t, 1> dimPcm{};
    pcmSpace.getSimpleExtentDims(dimPcm.data(), nullptr);
    std::vector<double> pcmVec(dimPcm[0]);
    pcm.read(pcmVec.data(), PredType::NATIVE_DOUBLE, pcmSpace, pcmSpace);

    // Get data for total cross-section
    DataSpace sigTotSpace = sigTot.getSpace();
    std::array<hsize_t, 1> dimSigTot{};
    sigTotSpace.getSimpleExtentDims(dimSigTot.data(), nullptr);
    std::vector<double> sigTotVec(dimSigTot[0]);
    sigTot.read(sigTotVec.data(), PredType::NATIVE_DOUBLE, sigTotSpace, sigTotSpace);

    // Get data for angular cross-section
    DataSpace sigSpace = sig.getSpace();
    std::array<hsize_t, 2> dimSig{};
    sigSpace.getSimpleExtentDims(dimSig.data(), nullptr);
    hsize_t dims = dimSig[0] * dimSig[1];
    std::vector<double> sigAngular(dims);
    sig.read(sigAngular.data(), PredType::NATIVE_DOUBLE, sigSpace, sigSpace);

    // Perform interpolation for angles
    nuchic::Interp2D interp(pcmVec, m_theta, sigAngular);
    interp.BicubicSpline();

    std::vector<double> theta(pcmVec.size()*m_cdf.size());
    constexpr double accuracy = 1E-6;
    for(size_t i = 0; i < pcmVec.size(); ++i) {
        for(size_t j = 0; j < m_cdf.size(); ++j) {
            auto func = [interp, pcmVec, this, i, j](double x){
                return interp(pcmVec[i], x) - this -> m_cdf[j];
            };
            nuchic::Brent brent(func, accuracy);
            if(j != m_cdf.size() - 1)
                try{
                    theta[i*m_cdf.size() + j] = brent.CalcRoot(m_theta.front(), m_theta.back());
                } catch (std::domain_error &e) {
                    theta[i*m_cdf.size() + j] = *m_theta.begin()/sigAngular[i*180 + j]*m_cdf[j];
                }
            else
                theta[i*m_cdf.size() + j] = *m_theta.end();
        }
    }

    if(samePID) {
        m_pcmPP = pcmVec;
        m_xsecPP = sigTotVec;
        m_crossSectionPP.SetData(pcmVec, sigTotVec);
        m_crossSectionPP.CubicSpline();
        m_thetaDistPP.SetData(pcmVec, m_cdf, theta);
        m_thetaDistPP.BicubicSpline();
    } else {
        m_pcmNP = pcmVec;
        m_xsecNP = sigTotVec;
        m_crossSectionNP.SetData(pcmVec, sigTotVec);
        m_crossSectionNP.CubicSpline();
        m_thetaDistNP.SetData(pcmVec, m_cdf, theta);
        m_thetaDistNP.BicubicSpline();
    }
}

double GeantInteractions::CrossSection(const Particle& particle1,
                                       const Particle& particle2) const {
    bool samePID = particle1.ID() == particle2.ID();
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM);
    // Generate outgoing momentum
    const double pcm = p1CM.Vec3().Magnitude();

    try {
        if(samePID) {
            return m_crossSectionPP(pcm/1_GeV);
        } else {
            return m_crossSectionNP(pcm/1_GeV);
        }
    } catch (std::domain_error &e) {
        spdlog::debug("Using Nasa Interaction");
        double s = (p1Lab+p2Lab).M2();
        double plab = sqrt(pow(s, 2)/(4*pow(Constant::mN, 2)) - s);
        return InteractionComponent::CrossSectionLab(samePID, plab);
    }
}

std::vector<Particle> GeantInteractions::GenerateFinalState(const Particle &particle1,
                                                            const Particle &particle2) const {

    // Boost to center of mass frame
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    // FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    // FourVector p1CM = p1Lab.Boost(-boostCM);

    std::vector<Particle> results = {particle1, particle2};
    std::array<double, 2> rans{};
    Random::Instance().Generate(rans, 0.0, 1.0);

    bool samePID = particle1.ID() == particle2.ID();
    double pcm = particle1.Momentum().Vec3().Magnitude();
    double theta = CrossSectionAngle(samePID, pcm/1_GeV, rans[0]);
    double phi = 2*M_PI*rans[1];

    auto mom = ThreeVector(ToCartesian({pcm, theta, phi}));
    auto p1Out = FourVector(mom[0], mom[1], mom[2], particle1.E());
    auto p2Out = FourVector(-mom[0], -mom[1], -mom[2], particle1.E());

    // Set momentum and boost back to lab frame
    results[0].SetMomentum(p1Out.Boost(boostCM));
    results[1].SetMomentum(p2Out.Boost(boostCM));

    return results;
}

double GeantInteractions::CrossSectionAngle(bool samePID, const double& energy,
                                            const double& ran) const {
    try{
        if(samePID) return m_thetaDistPP(energy, ran);
        else return m_thetaDistNP(energy, ran);
    } catch(std::domain_error &e) {
        spdlog::debug("Using flat angular distribution");
        return acos(2*ran-1);
    }
}

double NasaInteractions::CrossSection(const Particle& particle1,
                                       const Particle& particle2) const {
    bool samePID = particle1.ID() == particle2.ID();
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    // Generate outgoing momentum
    double s = (p1Lab+p2Lab).M2();
    double plab = sqrt(pow(s,2)/(4.0*pow(Constant::mN, 2))-s);
    return CrossSectionLab(samePID,plab); 
}

std::vector<Particle> NasaInteractions::GenerateFinalState(const Particle &particle1,
                                                           const Particle &particle2) const {
    // Boost to center of mass frame
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    // FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    // FourVector p1CM = p1Lab.Boost(-boostCM);

    std::vector<Particle> results = {particle1, particle2};
    std::array<double, 2> rans{};
    Random::Instance().Generate(rans, 0.0, 1.0);
    double ctheta = 2*rans[0]-1;
    double stheta = sqrt(1-ctheta*ctheta);
    double phi = 2*M_PI*rans[1];
   
    double pcm = particle1.Momentum().Vec3().Magnitude();
    auto p1Out = FourVector(pcm*stheta*cos(phi), pcm*stheta*sin(phi),
                            pcm*ctheta, particle1.E());
    auto p2Out = FourVector(-pcm*stheta*cos(phi), -pcm*stheta*sin(phi),
                            -pcm*ctheta, particle1.E());

    // Set momentum and boost back to lab frame
    results[0].SetMomentum(p1Out.Boost(boostCM));
    results[1].SetMomentum(p2Out.Boost(boostCM));

    return results;
}

std::vector<Particle> ConstantInteractions::GenerateFinalState(const Particle &particle1,
                                                               const Particle &particle2) const {
    // Boost to center of mass frame
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    // FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    // FourVector p1CM = p1Lab.Boost(-boostCM);

    std::vector<Particle> results = {particle1, particle2};
    std::array<double, 2> rans{};
    Random::Instance().Generate(rans, 0.0, 1.0);
    double ctheta = 2*rans[0]-1;
    double stheta = sqrt(1-ctheta*ctheta);
    double phi = 2*M_PI*rans[1];
   
    double pcm = particle1.Momentum().Vec3().Magnitude();
    auto p1Out = FourVector(pcm*stheta*cos(phi), pcm*stheta*sin(phi),
                            pcm*ctheta, particle1.E());
    auto p2Out = FourVector(-pcm*stheta*cos(phi), -pcm*stheta*sin(phi),
                            -pcm*ctheta, particle1.E());

    // Set momentum and boost back to lab frame
    results[0].SetMomentum(p1Out.Boost(boostCM));
    results[1].SetMomentum(p2Out.Boost(boostCM));

    return results;
}

PionNucleon::PionNucleon(const YAML::Node &config) {
    auto foldername = config["PionData"].as<std::string>();

    spdlog::info("Loading pion-nucleon cross sections from {}", foldername);
    m_cross_sections[{PID::pion0(), PID::neutron()}] = LoadData(foldername+filenames[0],
                     {{PID::pion0(), PID::neutron()}, {PID::pionm(), PID::proton()}});
    m_cross_sections[{PID::pion0(), PID::proton()}] = LoadData(foldername+filenames[1],
                     {{PID::pion0(), PID::proton()}, {PID::pionp(), PID::neutron()}});
    m_cross_sections[{PID::pionm(), PID::neutron()}] = LoadData(foldername+filenames[2],
                     {{PID::pionm(), PID::neutron()}});
    m_cross_sections[{PID::pionm(), PID::proton()}] = LoadData(foldername+filenames[3],
                     {{PID::pionm(), PID::proton()}, {PID::pion0(), PID::neutron()}});
    m_cross_sections[{PID::pionp(), PID::neutron()}] = LoadData(foldername+filenames[4],
                     {{PID::pionp(), PID::neutron()}, {PID::pion0(), PID::proton()}});
    m_cross_sections[{PID::pionp(), PID::proton()}] = LoadData(foldername+filenames[5],
                     {{PID::pionp(), PID::proton()}});
}

PionNucleon::cross_section PionNucleon::LoadData(const std::string &filename,
                                                 const std::vector<std::pair<PID, PID>> &pids) const {
    // Initialize the cross_section struct 
    cross_section results;
    results.m_pids = pids;

    // Initialize theta vector
    constexpr double thetaMin = 0.1_deg;
    constexpr double thetaMax = 179.9_deg;
    constexpr size_t nTheta = 180;
    results.m_theta = Linspace(thetaMin, thetaMax, nTheta);

    constexpr double cdfMin = -3;
    constexpr double cdfMax = 0;
    constexpr size_t nCDF = 200;
    results.m_cdf = Logspace(cdfMin, cdfMax, nCDF);

    // Load the data from file
    std::set<double> ecm, theta;
    std::vector<double> xsec1, xsec2;
    std::string line;
    std::ifstream data(filename); 
    if(!data.is_open()) {
        throw std::runtime_error(fmt::format("File Error: Could not open {}\n", filename));
    }
    std::getline(data, line); // Remove header line
    spdlog::debug("Reading pion data from: {}", filename);

    // Read in file
    while(std::getline(data, line)) {
        std::vector<std::string> tokens;
        tokenize(line, tokens);
        ecm.insert(std::stod(tokens[0]));
        theta.insert(std::stod(tokens[1])*1.0_deg);
        xsec1.emplace_back(std::stod(tokens[2]));
        if(tokens.size() == 4) {
            xsec2.emplace_back(std::stod(tokens[3]));
        } else {
            xsec2.emplace_back(0);
        }
    }
    data.close();
    spdlog::trace("Finished loading data for: {}", filename);

    // Calculate and invert the CDF numerically
    results.m_energies = {ecm.begin(), ecm.end()};
    results.m_angles = {theta.begin(), theta.end()};
    results.CalcCDF(xsec1, xsec2);

    return results;
}

double nuchic::PionNucleon::CrossSection(const Particle &pion, const Particle &nucleon) const {
    auto xsecs = CrossSections(pion, nucleon);
    return std::accumulate(xsecs.begin(), xsecs.end(), 0);
}

std::vector<double> nuchic::PionNucleon::CrossSections(const Particle &pion, const Particle &nucleon) const {
    double ecm = (pion.Momentum() + nucleon.Momentum()).M();
    spdlog::debug("ECM = {}, pids = ({}, {})", ecm, int(pion.ID()), int(nucleon.ID()));
    auto current_xsec = m_cross_sections.at({pion.ID(), nucleon.ID()});
    std::vector<double> results(current_xsec.m_cross_sections.size());
    for(size_t i = 0; i < results.size(); ++i) {
        try {
            results[i] = current_xsec.m_cross_sections[i](ecm);
        } catch (std::domain_error&) {
            results[i] = 0;
        }
    }
    return results;
}

std::vector<Particle> nuchic::PionNucleon::GenerateFinalState(const Particle &pion, const Particle &nucleon) const {
    // Boost to center of mass frame
    ThreeVector boostCM = (pion.Momentum() + nucleon.Momentum()).BoostVector();
    double ecm = (pion.Momentum() + nucleon.Momentum()).M();
    std::pair<PID, PID> init{pion.ID(), nucleon.ID()};
    size_t finalIdx = SelectFinalState(init, ecm, Random::Instance().Uniform(0.0, 1.0));
    auto pids = m_cross_sections.at(init).m_pids[finalIdx];

    std::vector<Particle> outgoing;
    std::array<double, 2> rans{};
    Random::Instance().Generate(rans, 0.0, 1.0);

    double pcm = pion.Momentum().Vec3().Magnitude();
    double theta = GenerateAngle(init, finalIdx, ecm, rans[0]);
    double phi = 2*M_PI*rans[1];

    auto mom = ThreeVector(ToCartesian({pcm, theta, phi}));
    auto pionOut = FourVector(mom[0], mom[1], mom[2], pion.E());
    auto nucleonOut = FourVector(-mom[0], -mom[1], -mom[2], nucleon.E());

    // Set momentum and boost back to lab frame
    outgoing.emplace_back(pids.first, pionOut.Boost(boostCM),
                          pion.Position());
    outgoing.emplace_back(pids.second, nucleonOut.Boost(boostCM),
                          nucleon.Position());

    return outgoing;
}

double nuchic::PionNucleon::GenerateAngle(const std::pair<PID, PID> &initial_state,
                                          const size_t &finalIdx,
                                          const double &energy,
                                          const double &ran) const {
    auto theta_dist = m_cross_sections.at(initial_state).m_theta_dist[finalIdx];
    try {
        return theta_dist(energy, ran);
    } catch (std::domain_error&) {
        return 0;
    }
}

size_t nuchic::PionNucleon::SelectFinalState(const std::pair<PID, PID>& initial_state,
                                             const double &energy,
                                             const double &ran) const {
    auto current_xsec = m_cross_sections.at(initial_state);
    std::vector<double> xsec(current_xsec.m_cross_sections.size());
    double total_xsec = 0;
    for(size_t i = 0; i < current_xsec.m_cross_sections.size(); ++i) {
        xsec[i] = current_xsec.m_cross_sections[i](energy);
        total_xsec += xsec[i];
    }

    size_t idx = 0;
    double prob = 0;
    for(const auto &ixsec : xsec) {
        prob += ixsec / total_xsec;
        if(prob > ran) break;
        ++idx;
    }
    return idx;
}

void nuchic::PionNucleon::cross_section::CalcCDF(const std::vector<double> &xsec1,
                                                 const std::vector<double> &xsec2) {
    const size_t nEnergy = m_energies.size();
    const size_t nAngle = m_angles.size();
    m_sigma1.resize(nEnergy);
    m_sigma2.resize(nEnergy);
    std::vector<double> theta1(nEnergy*m_cdf.size());
    std::vector<double> theta2(nEnergy*m_cdf.size());

    for(size_t iEnergy = 0; iEnergy < nEnergy; ++iEnergy) {
        m_sigma1[iEnergy] = 0;
        m_sigma2[iEnergy] = 0;

        // Integrate over angles to obtain total cross-section
        for(size_t iAngle = 0; iAngle < nAngle; ++iAngle) {
            m_sigma1[iEnergy] += xsec1[iEnergy*nAngle+iAngle]*sin(m_angles[iAngle])*10*2*M_PI;
            m_sigma2[iEnergy] += xsec2[iEnergy*nAngle+iAngle]*sin(m_angles[iAngle])*10*2*M_PI;
        }

        spdlog::trace("E = {}, sigma1 = {}, sigma2 = {}", m_energies[iEnergy], m_sigma1[iEnergy], m_sigma2[iEnergy]);

        // Calculate the angular cdf
        std::vector<double> cdf1(nAngle), cdf2(nAngle);
        for(size_t iAngle = 0; iAngle < nAngle; ++iAngle) {
            double sigma1 = xsec1[iEnergy*nAngle+iAngle]*sin(m_angles[iAngle])*10*2*M_PI;
            cdf1[iAngle] = sigma1/m_sigma1[iEnergy];
            if(iAngle != 0) cdf1[iAngle] += cdf1[iAngle-1];

            double sigma2 = xsec2[iEnergy*nAngle+iAngle]*sin(m_angles[iAngle])*10*2*M_PI;
            cdf2[iAngle] = sigma2/m_sigma2[iEnergy];
            if(iAngle != 0) cdf2[iAngle] += cdf2[iAngle-1];

            spdlog::trace("Angle = {}, cdf1 = {}, cdf2 = {}", m_angles[iAngle], cdf1[iAngle], cdf2[iAngle]);
        }

        // Create interpolation function
        nuchic::Interp1D interpCDF1(m_angles, cdf1);
        interpCDF1.CubicSpline();
        nuchic::Interp1D interpCDF2(m_angles, cdf2);
        interpCDF2.CubicSpline();

        // Invert the cdf for a fixed energy
        size_t idx = 0;
        constexpr double accuracy = 1E-6;
        for(const auto & cdf : m_cdf) {
            // Invert xsec1 cdf
            auto func1 = [interpCDF1, &cdf](double x) {
                return interpCDF1(x) - cdf;
            };
            nuchic::Brent brent1(func1, accuracy);
            try{
                theta1[iEnergy*m_cdf.size() + idx] = brent1.CalcRoot(m_theta.front(), m_theta.back());
            } catch (std::domain_error&) {
                theta1[iEnergy*m_cdf.size() + idx] = m_angles.front() / cdf1.front() * cdf;
            }

            // Invert xsec2 cdf
            auto func2 = [interpCDF2, &cdf](double x) {
                return interpCDF2(x) - cdf;
            };
            nuchic::Brent brent2(func2, accuracy);
            try{
                theta2[iEnergy*m_cdf.size() + idx] = brent2.CalcRoot(m_theta.front(), m_theta.back());
            } catch (std::domain_error&) {
                theta2[iEnergy*m_cdf.size() + idx] = m_angles.front() / cdf2.front() * cdf;
            }
            ++idx;
        }
    }

    m_cross_sections.emplace_back(m_energies, m_sigma1);
    m_cross_sections.back().CubicSpline();
    m_cross_sections.emplace_back(m_energies, m_sigma2);
    m_cross_sections.back().CubicSpline();
    m_theta_dist.emplace_back(m_energies, m_cdf, theta1);
    m_theta_dist.back().BicubicSpline();
    m_theta_dist.emplace_back(m_energies, m_cdf, theta2);
    m_theta_dist.back().BicubicSpline();
}
