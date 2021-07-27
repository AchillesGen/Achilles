#include "nuchic/LeptonicMapper.hh"
#include "nuchic/FourVector.hh"

using ATOOLS::Vec4D;
using nuchic::LeptonicMapper;
using nuchic::SherpaMapper;
using nuchic::FourVector;

void LeptonicMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double>&) const {
    // p_out = p_in - q_vec
    point[2] = point[1] - point[0];
}

double LeptonicMapper::GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const {
    return 1.0;
}

void SherpaMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    std::vector<Vec4D> mom(point.size());
    mom[0] = Vec4D(point[0][0], point[0][1], point[0][2], point[0][3]);
    mom[1] = Vec4D(point[1][0], point[1][1], point[1][2], point[1][3]);
    sherpa_mapper -> GeneratePoint(mom, rans);
    for(size_t i = 2; i < point.size(); ++i) {
        point[i] = FourVector(mom[i][0], mom[i][1], mom[i][2], mom[i][3]);
    }
}

double SherpaMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    std::vector<Vec4D> mom{};
    for(const auto &pt : point)
        mom.emplace_back(pt[0], pt[1], pt[2], pt[3]);
    return sherpa_mapper -> GenerateWeight(mom, rans);
}
