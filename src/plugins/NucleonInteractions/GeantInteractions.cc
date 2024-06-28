#include <iostream>
#include <map>

#include "spdlog/spdlog.h"

#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/GeantInteractions.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"

using namespace H5;
using namespace achilles;

REGISTER_INTERACTION(GeantInteractions);

GeantInteractions::GeantInteractions(const std::string &filename) {
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

void GeantInteractions::LoadData(bool samePID, const Group &group) {
    // Load datasets
    DataSet pcm(group.openDataSet("pcm"));
    DataSet sigTot(group.openDataSet("sigtot"));
    DataSet sig(group.openDataSet("sig"));

    // Get data for center of momentum
    DataSpace pcmSpace = pcm.getSpace();
    hsize_t dimPcm[1];
    pcmSpace.getSimpleExtentDims(dimPcm, nullptr);
    std::vector<double> pcmVec(dimPcm[0]);
    pcm.read(pcmVec.data(), PredType::NATIVE_DOUBLE, pcmSpace, pcmSpace);

    // Get data for total cross-section
    DataSpace sigTotSpace = sigTot.getSpace();
    hsize_t dimSigTot[1];
    sigTotSpace.getSimpleExtentDims(dimSigTot, nullptr);
    std::vector<double> sigTotVec(dimSigTot[0]);
    sigTot.read(sigTotVec.data(), PredType::NATIVE_DOUBLE, sigTotSpace, sigTotSpace);

    // Get data for angular cross-section
    DataSpace sigSpace = sig.getSpace();
    hsize_t dimSig[2];
    sigSpace.getSimpleExtentDims(dimSig, nullptr);
    hsize_t dims = dimSig[0] * dimSig[1];
    std::vector<double> sigAngular(dims);
    sig.read(sigAngular.data(), PredType::NATIVE_DOUBLE, sigSpace, sigSpace);

    // Perform interpolation for angles
    achilles::Interp2D interp;
    interp.BicubicSpline(pcmVec, m_theta, sigAngular);

    std::vector<double> theta(pcmVec.size() * m_cdf.size());
    constexpr double accuracy = 1E-6;
    for(size_t i = 0; i < pcmVec.size(); ++i) {
        for(size_t j = 0; j < m_cdf.size(); ++j) {
            auto func = [interp, pcmVec, this, i, j](double x) {
                return interp(pcmVec[i], x) - this->m_cdf[j];
            };
            achilles::Brent brent(func, accuracy);
            if(j != m_cdf.size() - 1) try {
                    theta[i * m_cdf.size() + j] = brent.CalcRoot(m_theta.front(), m_theta.back());
                } catch(std::domain_error &e) {
                    theta[i * m_cdf.size() + j] =
                        *m_theta.begin() / sigAngular[i * 180 + j] * m_cdf[j];
                }
            else
                theta[i * m_cdf.size() + j] = *m_theta.end();
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

double GeantInteractions::CrossSection(const Particle &particle1, const Particle &particle2) const {
    bool samePID = particle1.ID() == particle2.ID();
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM);
    // Generate outgoing momentum
    const double pcm = p1CM.Vec3().Magnitude();

    try {
        if(samePID) {
            return m_crossSectionPP(pcm / 1_GeV);
        } else {
            return m_crossSectionNP(pcm / 1_GeV);
        }
    } catch(std::domain_error &e) {
        spdlog::debug("Using Nasa Interaction");
        double s = (p1Lab + p2Lab).M2();
        double plab = sqrt(pow(s, 2) / (4 * pow(Constant::mN, 2)) - s);
        return Interactions::CrossSectionLab(samePID, plab);
    }
}

ThreeVector GeantInteractions::MakeMomentum(bool samePID, const double &pcm,
                                            const std::array<double, 2> &rans) const {
    double pR = pcm;
    double pTheta = CrossSectionAngle(samePID, pcm / 1_GeV, rans[0]);
    double pPhi = 2 * M_PI * rans[1];

    return ThreeVector(ToCartesian({pR, pTheta, pPhi}));
}

double GeantInteractions::CrossSectionAngle(bool samePID, const double &energy,
                                            const double &ran) const {
    try {
        if(samePID)
            return m_thetaDistPP(energy, ran);
        else
            return m_thetaDistNP(energy, ran);
    } catch(std::domain_error &e) {
        spdlog::debug("Using flat angular distribution");
        return acos(2 * ran - 1);
    }
}
