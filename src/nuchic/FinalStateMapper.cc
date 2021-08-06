#include "nuchic/FinalStateMapper.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"
#include "nuchic/Units.hh"

using ATOOLS::Vec4D;
using nuchic::TwoBodyMapper;
using nuchic::SherpaMapper;
using nuchic::FourVector;

void TwoBodyMapper::GeneratePoint(std::vector<FourVector> &mom, const std::vector<double> &rans) const {
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
    auto rotMat = mom0.AlignZ();
    auto cosT = dCos*rans[0] - 1;
    auto sinT = sqrt(1 - cosT*cosT);
    auto phi = dPhi*rans[1];
    auto E1 = sqrts/2*(1 + s2/s - s3/s);
    auto E2 = sqrts/2*(1 + s3/s - s2/s);
    auto beta = sqrt(1- 2*(s2+s3)/s + pow(s2 - s3, 2)/s/s);
    auto pCM = sqrts/2*beta;

    mom[2] = {E1, pCM*sinT*cos(phi), pCM*sinT*sin(phi), pCM*cosT};
    mom[3] = {E2, -pCM*sinT*cos(phi), -pCM*sinT*sin(phi), -pCM*cosT};

    mom[2] = mom[2].RotateBack(rotMat).Boost(boostVec);
    mom[3] = mom[3].RotateBack(rotMat).Boost(boostVec);

    Mapper<nuchic::FourVector>::Print(__PRETTY_FUNCTION__, mom, rans);
}

double TwoBodyMapper::GenerateWeight(const std::vector<FourVector> &mom, std::vector<double> &rans) const {
    auto boostVec = (mom[0] + mom[1]).BoostVector();
    auto mom0 = mom[0].Boost(-boostVec);
    auto rotMat = mom0.AlignZ();
    auto p2 = mom[2].Boost(-boostVec).Rotate(rotMat);
    rans[0] = (p2.CosTheta() + 1)/dCos;
    rans[1] = p2.Phi()/dPhi;

    auto pcm = mom0.P();
    auto ecm = (mom[0] + mom[1]).E();

    auto factor = pcm/ecm*mom[3].E()/mom[1].E();
    Mapper<nuchic::FourVector>::Print(__PRETTY_FUNCTION__, mom, rans);
    spdlog::trace("  Factor: {}", factor);

    return 1.0/dCos/dPhi/factor;
}

void SherpaMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    std::vector<Vec4D> mom(point.size());
    mom[0] = Vec4D(point[0][0]/1_GeV, point[0][1]/1_GeV, point[0][2]/1_GeV, point[0][3]/1_GeV);
    mom[1] = Vec4D(point[1][0]/1_GeV, point[1][1]/1_GeV, point[1][2]/1_GeV, point[1][3]/1_GeV);
    sherpa_mapper -> GeneratePoint(mom, rans);
    for(size_t i = 2; i < point.size(); ++i) {
        point[i] = FourVector(mom[i][0]*1_GeV, mom[i][1]*1_GeV, mom[i][2]*1_GeV, mom[i][3]*1_GeV);
    }
}

double SherpaMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    std::vector<Vec4D> mom{};
    for(const auto &pt : point)
        mom.emplace_back(pt[0]/1_GeV, pt[1]/1_GeV, pt[2]/1_GeV, pt[3]/1_GeV);
    return sherpa_mapper -> GenerateWeight(mom, rans);
}
