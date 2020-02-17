#include <iostream>
#include <map>

#include "nuchic/FourVector.hh"
#include "nuchic/Interactions.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"
#include "nuchic/Utilities.hh"

#define FM 1.0
#define TO_NB 1e6
#define MEV 1.0
#define GEV 1.e3
#define C 2.9979245858E8*M
#define HBARC 197.3269788*MEV*FM
#define MP 938
#define MN 939
#define HBARC2 3.8937966E8

using namespace H5;

const std::map<std::string, double> HZETRN = {
    {"a", 5.0 * MEV},
    {"b", 0.199/sqrt(MEV)},
    {"c", 0.451 * pow(MEV, -0.258)},
    {"d", 25.0 * MEV},
    {"e", 134.0 * MEV},
    {"f", 1.187 * pow(MEV, -0.35)},
    {"g", 0.1 * MEV},
    {"h", 0.282 * MEV}
};

const std::map<std::string, double> PDG = {
    {"Zpp", 33.45},  //mb
    {"Zpn", 35.80},  //mb
    {"Y1pp", 42.53}, //mb
    {"Y1pn", 40.15}, //mb
    {"Y2pp", 33.34}, //mbs
    {"Y2pn", 30.00}, //mb
    {"B", 0.308},    //mb
    {"s1", 1.0 * pow(GEV, 2)},
    {"s0", pow(5.38 * GEV, 2)},
    {"n1", 0.458},
    {"n2", 0.545}
};

const std::map<std::string, double> JWN = {
    {"gamma", 52.5 * pow(GEV, 0.16)}, //mb
    {"alpha", 0.00369 / MEV},
    {"beta", 0.00895741 * pow(MEV, -0.8)}
};

double nuchic::CrossSection(bool samePID, const double& pcm) {
    // TODO: Implement GEANT4 cross-section or something else
    throw std::domain_error("Invalid energy!");
}

double nuchic::CrossSectionLab(bool samePID, const double& pLab) noexcept {
    const double tLab = sqrt(pow(pLab, 2) + pow(MN, 2)) - MN;
    if(samePID) {
        if(pLab < 1.8 * GEV) {
            if(tLab >= 25 * MEV)
                return (1.0+HZETRN.at("a")/tLab) * (40+109.0*std::cos(HZETRN.at("b")*sqrt(tLab))
                        * exp(-HZETRN.at("c")*pow(tLab-HZETRN.at("d"), 0.258)));
            else
                return exp(6.51*exp(-pow(tLab/HZETRN.at("e"), 0.7)));
        } else if(pLab <= 4.7 * GEV) {
            return JWN.at("gamma")/pow(pLab, 0.16);
        } else {
            double ecm2 = 2*MN*(MN+sqrt(pow(pLab, 2) + pow(MN, 2)));
            return PDG.at("Zpp") + PDG.at("B")*pow(log(ecm2/PDG.at("s0")), 2)
                + PDG.at("Y1pp")*pow(PDG.at("s1")/ecm2, PDG.at("n1"))
                - PDG.at("Y2pp")*pow(PDG.at("s1")/ecm2, PDG.at("n2"));
        }
    } else {
        if(pLab < 0.5 * 1e3) {
            if(tLab >= 0.1 * MEV)
                return 38.0 + 12500.0*exp(-HZETRN.at("f")*pow(tLab-HZETRN.at("g"), 0.35));
            else
                return 26000 * exp(-pow(tLab/HZETRN.at("h"), 0.3));
        } else if(pLab <= 2.0 * 1e3) {
            return 40 + 10*cos(JWN.at("alpha")*pLab - 0.943)
                * exp(-JWN.at("beta")*pow(pLab, 0.8)+2);
        } else {
            double ecm2 = 2*MN*(MN+sqrt(pow(pLab, 2) + pow(MN, 2)));
            return PDG.at("Zpn") + PDG.at("B")*pow(log(ecm2/PDG.at("s0")), 2)
                + PDG.at("Y1pn")*pow(PDG.at("s1")/ecm2, PDG.at("n1"))
                - PDG.at("Y2pn")*pow(PDG.at("s1")/ecm2, PDG.at("n2"));
        }
    }
}

double nuchic::CrossSectionAngle(bool samePID, const double& pcm, const double& ran) {
    // For testing right now
    // TODO: Implement GEANT4 angular distribution or something else
    return std::acos(2*ran - 1);
}

nuchic::ThreeVector nuchic::MakeMomentumAngular(bool samePID, const double& p1CM, const double& pcm,
        const std::array<double, 2>& rans) {
    double pR = p1CM;
    double pTheta = nuchic::CrossSectionAngle(samePID, pcm, rans[0]);
    double pPhi = 2*M_PI*rans[1];

    return nuchic::ThreeVector(nuchic::ToCartesian({pR, pTheta, pPhi}));
}

double nuchic::Interactions::CrossSectionLab(bool samePID, const double& pLab) const noexcept {
    const double tLab = sqrt(pow(pLab, 2) + pow(MN, 2)) - MN;
    if(samePID) {
        if(pLab < 1.8 * GEV) {
            if(tLab >= 25 * MEV)
                return (1.0+HZETRN.at("a")/tLab) * (40+109.0*std::cos(HZETRN.at("b")*sqrt(tLab))
                        * exp(-HZETRN.at("c")*pow(tLab-HZETRN.at("d"), 0.258)));
            else
                return exp(6.51*exp(-pow(tLab/HZETRN.at("e"), 0.7)));
        } else if(pLab <= 4.7 * GEV) {
            return JWN.at("gamma")/pow(pLab, 0.16);
        } else {
            double ecm2 = 2*MN*(MN+sqrt(pow(pLab, 2) + pow(MN, 2)));
            return PDG.at("Zpp") + PDG.at("B")*pow(log(ecm2/PDG.at("s0")), 2)
                + PDG.at("Y1pp")*pow(PDG.at("s1")/ecm2, PDG.at("n1"))
                - PDG.at("Y2pp")*pow(PDG.at("s1")/ecm2, PDG.at("n2"));
        }
    } else {
        if(pLab < 0.5 * GEV) {
            if(tLab >= 0.1 * MEV)
                return 38.0 + 12500.0*exp(-HZETRN.at("f")*pow(tLab-HZETRN.at("g"), 0.35));
            else
                return 26000 * exp(-pow(tLab/HZETRN.at("h"), 0.3));
        } else if(pLab <= 2.0 * GEV) {
            return 40 + 10*cos(JWN.at("alpha")*pLab - 0.943)
                * exp(-JWN.at("beta")*pow(pLab, 0.8)+2);
        } else {
            double ecm2 = 2*MN*(MN+sqrt(pow(pLab, 2) + pow(MN, 2)));
            return PDG.at("Zpn") + PDG.at("B")*pow(log(ecm2/PDG.at("s0")), 2)
                + PDG.at("Y1pn")*pow(PDG.at("s1")/ecm2, PDG.at("n1"))
                - PDG.at("Y2pn")*pow(PDG.at("s1")/ecm2, PDG.at("n2"));
        }
    }
}

REGISTER_INTERACTION(nuchic::GeantInteractions);

nuchic::GeantInteractions::GeantInteractions(const std::string& filename) : nuchic::Interactions(filename) {
    // Initialize theta vector
    m_theta = nuchic::Linspace(0.5, 179.5, 180);
    m_cdf = nuchic::Logspace(-3, 0, 200);

    // Read in the Geant4 hdf5 file and get the np and pp groups
    std::cout << "Loading Geant4 data..." << std::endl;
    H5File file(filename, H5F_ACC_RDONLY);
    Group dataNP(file.openGroup("np")); 
    Group dataPP(file.openGroup("pp"));

    // Get the datasets for np and load into local variables
    LoadData(false, dataNP);

    // Get the datasets for pp and load into local variables
    LoadData(true, dataPP);
}

void nuchic::GeantInteractions::LoadData(bool samePID, const Group& group) {
    // Load datasets
    DataSet pcm(group.openDataSet("pcm"));
    DataSet sigTot(group.openDataSet("sigtot"));
    DataSet sig(group.openDataSet("sig"));
  
    // Get data for center of momentum
    DataSpace pcmSpace = pcm.getSpace();
    hsize_t dimPcm[1];
    pcmSpace.getSimpleExtentDims(dimPcm, NULL);
    std::vector<double> pcmVec(dimPcm[0]);
    pcm.read(pcmVec.data(), PredType::NATIVE_DOUBLE, pcmSpace, pcmSpace);

    // Get data for total cross-section
    DataSpace sigTotSpace = sigTot.getSpace();
    hsize_t dimSigTot[1];
    sigTotSpace.getSimpleExtentDims(dimSigTot, NULL);
    std::vector<double> sigTotVec(dimSigTot[0]);
    sigTot.read(sigTotVec.data(), PredType::NATIVE_DOUBLE, sigTotSpace, sigTotSpace);

    // Get data for angular cross-section
    DataSpace sigSpace = sig.getSpace();
    hsize_t dimSig[2];
    sigSpace.getSimpleExtentDims(dimSig, NULL);
    hsize_t dims = dimSig[0] * dimSig[1];
    std::vector<double> sigAngular(dims);
    sig.read(sigAngular.data(), PredType::NATIVE_DOUBLE, sigSpace, sigSpace);

    // Perform interpolation for angles
    nuchic::Interp2D interp;
    interp.BicubicSpline(pcmVec, m_theta, sigAngular);

    std::vector<double> theta(pcmVec.size()*m_cdf.size());
    for(size_t i = 0; i < pcmVec.size(); ++i) {
        for(size_t j = 0; j < m_cdf.size(); ++j) {
            auto func = [interp, pcmVec, this, i, j](double x){
                return interp(pcmVec[i], x) - this -> m_cdf[j];
            };
            nuchic::Brent brent(func, 1e-6);
            if(j != m_cdf.size() - 1)
                try{
                    theta[i*m_cdf.size() + j] = brent.CalcRoot(0.5, 179.5);
                } catch (std::runtime_error &e) {
                    theta[i*m_cdf.size() + j] = 0.5/sigAngular[i*180 + j] * m_cdf[j];
                }
            else
                theta[i*m_cdf.size() + j] = 179.5;
        }
    }

    if(samePID) {
        m_pcmPP = pcmVec;
        m_xsecPP = sigTotVec;
        m_crossSectionPP.CubicSpline(pcmVec, sigTotVec);
        m_thetaDistPP.BicubicSpline(pcmVec, m_cdf, theta);
    } else {
        m_pcmNP = pcmVec;
        m_xsecNP = sigTotVec;
        m_crossSectionNP.CubicSpline(pcmVec, sigTotVec);
        m_thetaDistNP.BicubicSpline(pcmVec, m_cdf, theta);
    }
}

double nuchic::GeantInteractions::CrossSection(const Particle& particle1,
                                               const Particle& particle2) const {
    bool samePID = particle1.PID() == particle2.PID();
    const nuchic::FourVector totalMomentum = particle1.Momentum() + particle2.Momentum();
    const double pcm = particle1.Momentum().Vec3().Magnitude() * MN / totalMomentum.M();


    try {
        if(samePID)
            return m_crossSectionPP(pcm/GEV);
        else
            return m_crossSectionNP(pcm/GEV);
    } catch (std::domain_error &e) {
        nuchic::FourVector momentum1 = particle1.Momentum();
        nuchic::FourVector momentum2 = particle2.Momentum();
        nuchic::FourVector labMomentum = momentum1.Boost(-momentum2.BoostVector());
        return nuchic::CrossSectionLab(samePID, labMomentum.Vec3().Magnitude()/GEV);
    }
}

nuchic::ThreeVector nuchic::GeantInteractions::MakeMomentum(bool samePID, const double& p1CM,
                                                            const double& pcm,
                                                            const std::array<double, 2>& rans) const {
    double pR = p1CM;
    double pTheta = CrossSectionAngle(samePID, pcm/GEV, rans[0]);
    double pPhi = 2*M_PI*rans[1];

    return nuchic::ThreeVector(nuchic::ToCartesian({pR, pTheta, pPhi}));
}

double nuchic::GeantInteractions::CrossSectionAngle(bool samePID, const double& energy,
        const double& ran) const {
    if(samePID) return m_thetaDistPP(energy, ran);
    else return m_thetaDistNP(energy, ran);
}
