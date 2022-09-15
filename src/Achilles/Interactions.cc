#include <fstream>
#include <map>

#include "Achilles/Potential.hh"
#include "spdlog/spdlog.h"

#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"
#include "Achilles/Random.hh"

using namespace achilles;
using namespace H5;

REGISTER_INTERACTION(GeantInteractions);
// REGISTER_INTERACTION(GeantInteractionsDt);
REGISTER_INTERACTION(NasaInteractions);
REGISTER_INTERACTION(ConstantInteractions);

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

double Interactions::CrossSectionLab(bool samePID, const double& pLab) const noexcept {
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

achilles::Interactions::MomentumPair Interactions::FinalizeMomentum(const Particle &particle1,
                                                                    const Particle &particle2,
                                                                    Potential*) const {

//    if(pot -> IsRelativistic())
//        throw std::runtime_error(fmt::format("{} is not compatible with a relativistic potential.",
//                                             Name()));
   
    // Boost to center of mass
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM);

    // Generate outgoing momentum
    bool samePID = particle1.ID() == particle2.ID();
    const double pcm = p1CM.Vec3().Magnitude();
    std::array<double, 2> rans{};
    Random::Instance().Generate(rans, 0.0, 1.0);
    ThreeVector momentum = MakeMomentum(samePID, pcm, rans);

    FourVector p1Out = FourVector(p1CM.E(), momentum[0], momentum[1], momentum[2]);
    FourVector p2Out = FourVector(p1CM.E(), -momentum[0], -momentum[1], -momentum[2]);

    // Boost back to lab frame
    p1Out = p1Out.Boost(boostCM);
    p2Out = p2Out.Boost(boostCM);

    return {p1Out, p2Out};
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

    // Clean-up
    dataNP.close();
    dataPP.close();
    file.close();
    H5Library::close();
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
    achilles::Interp2D interp(pcmVec, m_theta, sigAngular);
    interp.BicubicSpline();

    std::vector<double> theta(pcmVec.size()*m_cdf.size());
    constexpr double accuracy = 1E-6;
    for(size_t i = 0; i < pcmVec.size(); ++i) {
        for(size_t j = 0; j < m_cdf.size(); ++j) {
            auto func = [interp, pcmVec, this, i, j](double x){
                return interp(pcmVec[i], x) - this -> m_cdf[j];
            };
            achilles::Brent brent(func, accuracy);
            if(j != m_cdf.size() - 1)
                try{
                    theta[i*m_cdf.size() + j] = brent.CalcRoot(m_theta.front(), m_theta.back());
                } catch (std::domain_error &e) {
                    theta[i*m_cdf.size() + j] = m_theta.front()/sigAngular[i*180 + j]*m_cdf[j];
                }
            else
                theta[i*m_cdf.size() + j] = m_theta.back();
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

    // Clean-up
    pcmSpace.close();
    sigTotSpace.close();
    sigSpace.close();
    pcm.close();
    sigTot.close();
    sig.close();
}

double GeantInteractions::CrossSection(const Particle& particle1,
                                       const Particle& particle2) const {
    bool samePID = particle1.ID() == particle2.ID();
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum(); //, p2Lab = particle2.Momentum();
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
        spdlog::trace("Using Nasa Interaction");
        // double s = (p1Lab+p2Lab).M2();
        double s = (particle1.Momentum()+particle2.Momentum()).M2();
        double smin = pow(particle1.Mass(), 2) + pow(particle2.Mass(), 2);
        double plab = sqrt(pow(s, 2)/smin - s);
        return Interactions::CrossSectionLab(samePID, plab);
    }
}

ThreeVector GeantInteractions::MakeMomentum(bool samePID,
                                            const double& pcm,
                                            const std::array<double, 2>& rans) const {
    double pR = pcm;
    double pTheta = CrossSectionAngle(samePID, pcm/1_GeV, rans[0]);
    double pPhi = 2*M_PI*rans[1];

    return ThreeVector(ToCartesian({pR, pTheta, pPhi}));
}

double GeantInteractions::CrossSectionAngle(bool samePID, const double& energy,
                                            const double& ran) const {
    try{
        if(samePID) return m_thetaDistPP(energy, ran);
        else return m_thetaDistNP(energy, ran);
    } catch(std::domain_error &e) {
        spdlog::trace("Using flat angular distribution");
        return acos(2*ran-1);
    }
}

/*
GeantInteractionsDt::GeantInteractionsDt(const YAML::Node& node) {
    auto filename_pp = node["GeantDataPP"].as<std::string>();
    auto filename_np = node["GeantDataNP"].as<std::string>();

    // Read in the Geant4 dsigma/dt files
    spdlog::info("GeantInteractionsDt: Loading Geant4 dsigma/dt_pp data from {0}.", filename_pp);
    LoadData(true, filename_pp);

    spdlog::info("GeantInteractionsDt: Loading Geant4 dsigma/dt_np data from {0}.", filename_np);
    LoadData(false, filename_np);
    throw;
}

void GeantInteractionsDt::LoadData(bool samePID, const std::string &filename) {
    std::ifstream data(filename);

    // Read in CDFs
    const size_t nlines = samePID ? 40 : 39;
    auto cdf = ReadBlock(data, nlines);

    // Read in tmin and tmax
    auto tmin = ReadBlock(data, 2);
    auto tmax = ReadBlock(data, 1);

    // Read in pcm
    auto ecm2 = ReadBlock(data, 2);
    for(auto & e : ecm2)
        e = 4*(e*e + Constant::mN*Constant::mN);

    // Read in lab energy and max sigma
    auto elab = ReadBlock(data, 2);
    auto max_sig = ReadBlock(data, 2);

    // Read in total cross-section
    auto xsec = ReadBlock(data, 2);

    data.close();

    // Obtain the PDF from the CDF
    const size_t npts = cdf.size()/tmin.size();
    std::vector<double> pdf(ecm2.size()*npts);
    for(size_t i = 0; i < ecm2.size(); ++i) {
        double tmp = 0;
        for(size_t j = 0; j < npts; ++j) {
            if(j == 0) pdf[i*npts] = cdf[i*npts];
            else pdf[i*npts + j] = cdf[i*npts+j] - cdf[i*npts+j-1];

            tmp += pdf[i*npts + j];
        }
    }

    // Store data in class
    auto trange = Linspace(0, 1, nlines);
    std::reverse(trange.begin(), trange.end());
    if(samePID) {
        m_crossSectionPP.SetData(ecm2, xsec);
        m_crossSectionPP.CubicSpline();
        m_pdfPP.SetData(ecm2, trange, pdf);
        m_pdfPP.BicubicSpline();
    } else {
        m_crossSectionNP.SetData(ecm2, xsec);
        m_crossSectionNP.CubicSpline();
        m_pdfNP.SetData(ecm2, trange, pdf);
        m_pdfNP.BicubicSpline();
    }
}

std::vector<double> GeantInteractionsDt::ReadBlock(std::ifstream &file, size_t nlines) const {
    std::string line{};
    std::vector<std::string> tokens;
    for(size_t i = 0; i < nlines; ++i) {
        std::getline(file, line);
        tokenize(line, tokens);
    }
    std::vector<double> results;
    for(const auto& token : tokens)
        results.emplace_back(std::stod(token));

    return results;
}

double GeantInteractionsDt::CrossSection(const Particle &particle1,
                                         const Particle &particle2) const {
    bool samePID = particle1.ID() == particle2.ID(); 
    const double ecm2 = (particle1.Momentum() + particle2.Momentum()).M2()/1_GeV/1_GeV;

    try{
        if(samePID) {
            return m_crossSectionPP(ecm2);
        } else {
            return m_crossSectionNP(ecm2);
        }
    } catch(std::domain_error &e) {
        spdlog::warn("ecm2 {} not in table, using 0 for xsec", ecm2);
        return 0;
    }
}

// TODO: Finish implementing this function
achilles::Interactions::MomentumPair GeantInteractionsDt::FinalizeMomentum(const Particle &particle1,
                                                                         const Particle &particle2,
                                                                         std::shared_ptr<Potential> pot) const {
    // Generate outgoing momentum
    // bool samePID = particle1.ID() == particle2.ID();
    // const double ecm2 = (particle1.Momentum() + particle2.Momentum()).M2();
    // const double tmin = 2*particle1.Mass() + 2*particle2.Mass() - ecm2;
    auto p1 = particle1.Momentum(); 
    auto p2 = particle2.Momentum();
    // auto q_free = p1 + p2;
    p1.E() = pot->Hamiltonian(p1.P(), particle1.Radius());
    p2.E() = pot->Hamiltonian(p2.P(), particle2.Radius());
    // auto q = p1+p2;
    
    // Rotate so (p1 + p2) is along the z-axis
    // auto rotation  = q.AlignZ();
    // q = q.Rotate(rotation);

    // Generate random variables
    // std::array<double, 2> rans{};
    // Random::Instance().Generate(rans, 0.0, 1.0);
    // const double phi = 2*M_PI*rans[0];
    // const double t =  samePID ?
    //     (m_pdfPP(ecm2, rans[1])-1)*tmin : (m_pdfNP(ecm2, rans[1])-1)*tmin;

    // ThreeVector momentum = MakeMomentum(samePID, ecm2, rans);

    return {p1, p2};
}

// TODO: Finish implementing this function
//ThreeVector GeantInteractionsDt::MakeMomentum(bool samePID,
//                                              const double &pcm,
//                                              const std::array<double, 2> &rans) const {
//    double pR = pcm;
//}
*/

double NasaInteractions::CrossSection(const Particle& particle1,
                                      const Particle& particle2) const {
    bool samePID = particle1.ID() == particle2.ID();
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    // Generate outgoing momentum
    double s = (p1Lab+p2Lab).M2();
    double smin = pow(particle1.Mass() + particle2.Mass(), 2);
    double plab = sqrt(pow(s,2)/smin-s);
    return CrossSectionLab(samePID,plab); 
}

ThreeVector NasaInteractions::MakeMomentum(bool, const double& pcm,
                                           const std::array<double, 2>& rans) const {
    double pR = pcm;
    double pTheta = acos(2*rans[0]-1);
    double pPhi = 2*M_PI*rans[1];

    return ThreeVector(ToCartesian({pR, pTheta, pPhi}));
}
