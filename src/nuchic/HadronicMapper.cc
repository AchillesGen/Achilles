#include "nuchic/HadronicMapper.hh"
#include "nuchic/Constants.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"
#include "spdlog/spdlog.h"

using nuchic::QESpectralMapper;

void QESpectralMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    // Generate inital nucleon state
    const double mom = dp*rans[0];
    const double cosT = dCos*rans[1] - 1;
    const double sinT = sqrt(1 - cosT*cosT);
    const double phi = dPhi*rans[2];
    ThreeVector pmom = {mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT};

    const double det = pow(point[1].E(), 2) + mom*mom + 2*pmom*point[1].Vec3() + pow(Constant::mN, 2);
    const double emax = Constant::mN + point[1].E() - sqrt(det);
    const double energy = emax*rans[3] - 1e-8;

    point[HadronIdx()] = {Constant::mN - energy, mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT};
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
    spdlog::trace("  mom = {}", mom);
    spdlog::trace("  emax = {}", emax);
    spdlog::trace("  emax2 = {}", - point[1].E() - sqrt(det));
    spdlog::trace("  s = {}", (point[0] + point[1]).M2());
}

double QESpectralMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    rans[0] = point[HadronIdx()].P()/dp;
    rans[1] = (point[HadronIdx()].CosTheta()+1)/dCos;
    rans[2] = point[HadronIdx()].Phi()/dPhi;

    const double det = pow(point[1].E(), 2) + point[0].P2() 
                     + 2*point[0].Vec3()*point[1].Vec3() + pow(Constant::mN, 2);
    const double emax = Constant::mN + point[1].E() - sqrt(det);
    const double energy = Constant::mN - point[HadronIdx()].E();
    rans[3] = (energy + 1e-8)/emax;
    const double dE = emax;
    double wgt = 1.0/point[0].P2()/dp/dCos/dPhi/dE;
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
    spdlog::trace("  Weight: {}", wgt);

    return wgt;
}
