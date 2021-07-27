#include "nuchic/HadronicMapper.hh"
#include "nuchic/Constants.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"

using nuchic::QESpectralMapper;

void QESpectralMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    // Generate phase space
    double cosT = dCos*rans[0] - 1;
    double sinT = sqrt(1 - cosT * cosT);
    double phi = dPhi*rans[1];
    double p = dp*rans[2];

    ThreeVector tmp = { p*sinT*cos(phi), p*sinT*sin(phi), p*cosT };
    tmp += point[0].Vec3();
    
    double Epp = sqrt(pow(Constant::mN, 2) + tmp.P2());
    double Ep = Constant::mN + point[0].E() - Epp;
    point[point.size()-2] = {Ep, p*sinT*cos(phi), p*sinT*sin(phi), p*cosT};
    point[point.size()-1] = {tmp, Epp};
}

double QESpectralMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    rans[0] = point[point.size()-2].P()/dp;
    rans[1] = point[point.size()-2].CosTheta()/dCos;
    rans[2] = point[point.size()-2].Phi()/dPhi;

    // Pauli Blocking
    // TODO: Is this correct? This would assume a global Fermi gas,
    //       but this is then inconsistent with a local Fermi gas for the cascade
    // wt *= point.back().P() > 225 ? 1 : 0;

    return point[point.size()-2].E() > 0 ? dp*point[point.size()-2].P2()*dCos*dPhi : 0;
}
