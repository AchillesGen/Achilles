#include "Achilles/FinalStateMapper.hh"
#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Poincare.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Units.hh"
#include "Achilles/Utilities.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
using ATOOLS::Vec4D;
using achilles::SherpaMapper;
#endif
using achilles::TwoBodyMapper;
using achilles::FourVector;

void TwoBodyMapper::GeneratePoint(std::vector<FourVector> &mom, const std::vector<double> &rans) {
    // The momentum are given in the following order:
    // 1. Momentum of the initial hadron
    // 2. Momentum of the initial lepton
    // 3. Momentum of all outgoing parts of the leptonic tensor
    // 4. Momentum of all outgoing hadrons
    auto p01 = (mom[0] + mom[1]);
    auto s = p01.M2();
    auto sqrts = sqrt(s);
    auto boostVec = p01.BoostVector();
    auto mom0 = mom[0].Boost(-boostVec);
    Poincare zax(mom0, FourVector(1.,0.,0.,1.));
    auto cosT = dCos*rans[0] - 1;
    auto sinT = sqrt(1 - cosT*cosT);
    auto phi = dPhi*rans[1];
    auto E1 = sqrts/2*(1 + s2/s - s3/s);
    auto E2 = sqrts/2*(1 + s3/s - s2/s);
    auto lambda = sqrt(pow(s-s2-s3, 2) - 4*s2*s3);
    auto pCM = lambda/(2*sqrts);

    mom[2] = {E1, pCM*sinT*cos(phi), pCM*sinT*sin(phi), pCM*cosT};
    mom[3] = {E2, -pCM*sinT*cos(phi), -pCM*sinT*sin(phi), -pCM*cosT};

    zax.RotateBack(mom[2]);
    zax.RotateBack(mom[3]);

    mom[2] = mom[2].Boost(boostVec);
    mom[3] = mom[3].Boost(boostVec);

    Mapper<achilles::FourVector>::Print(__PRETTY_FUNCTION__, mom, rans);
    spdlog::trace("  MassCheck: {}", CheckMasses({mom[2], mom[3]}, {s2, s3}));
    spdlog::trace("  s = {}, lambda = {}", s, lambda);
}

double TwoBodyMapper::GenerateWeight(const std::vector<FourVector> &mom, std::vector<double> &rans) {
    auto boostVec = (mom[0] + mom[1]).BoostVector();
    auto mom0 = mom[0].Boost(-boostVec);
    auto rotMat = mom0.AlignZ();
    auto p2 = mom[2].Boost(-boostVec).Rotate(rotMat);
    rans[0] = (p2.CosTheta() + 1)/dCos;
    rans[1] = p2.Phi()/dPhi;

    auto pcm = p2.P();
    auto ecm = (mom[0] + mom[1]).M();

    auto factor = pcm/ecm/(16*M_PI*M_PI);
    auto wgt = 1.0/dCos/dPhi/factor;
    Mapper<achilles::FourVector>::Print(__PRETTY_FUNCTION__, mom, rans);
    spdlog::trace("  ct: {}", p2.CosTheta());
    spdlog::trace("  pcm: {}", pcm);
    spdlog::trace("  ecm: {}", ecm);
    spdlog::trace("  Weight: {}", wgt);

    return wgt;
}

#ifdef ACHILLES_SHERPA_INTERFACE
void SherpaMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) {
    std::vector<Vec4D> mom(point.size());
    mom[0] = Vec4D(point[0][0]/1_GeV, point[0][1]/1_GeV, point[0][2]/1_GeV, point[0][3]/1_GeV);
    mom[1] = Vec4D(point[1][0]/1_GeV, point[1][1]/1_GeV, point[1][2]/1_GeV, point[1][3]/1_GeV);
    sherpa_mapper -> GeneratePoint(mom, rans);
    for(size_t i = 2; i < point.size(); ++i) {
        point[i] = FourVector(mom[i][0]*1_GeV, mom[i][1]*1_GeV, mom[i][2]*1_GeV, mom[i][3]*1_GeV);
    }
}

double SherpaMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) {
    std::vector<Vec4D> mom{};
    for(const auto &pt : point)
        mom.emplace_back(pt[0]/1_GeV, pt[1]/1_GeV, pt[2]/1_GeV, pt[3]/1_GeV);
    return sherpa_mapper -> GenerateWeight(mom, rans);
}
#endif
