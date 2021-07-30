#include "nuchic/HadronicMapper.hh"
#include "nuchic/Constants.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"

using nuchic::QESpectralMapper;

void QESpectralMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    // Generate inital nucleon state
    const double mom = dp*rans[0];
    const double cosT = dCos*rans[1] - 1;
    const double sinT = sqrt(1 - cosT*cosT);
    const double phi = dPhi*rans[2];
    const double energy = dE*rans[3];  

    point[HadronIdx()] = {Constant::mN - energy, mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT};
}

double QESpectralMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    rans[0] = point[HadronIdx()].P()/dp;
    rans[1] = (point[HadronIdx()].CosTheta()+1)/dCos;
    rans[2] = point[HadronIdx()].Phi()/dPhi;
    rans[3] = -(point[HadronIdx()].E() - Constant::mN)/dE;

    return 1.0/point[HadronIdx()].P2()/dp/dCos/dPhi/dE;
}
