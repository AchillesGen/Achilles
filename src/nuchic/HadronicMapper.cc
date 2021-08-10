#include "nuchic/HadronicMapper.hh"
#include "nuchic/Constants.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"
#include "spdlog/spdlog.h"

using nuchic::QESpectralMapper;

void QESpectralMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    // Generate inital nucleon state
    constexpr double dp = 800;
    const double mom = dp*rans[0];
    // constexpr double lambda = 200;
    // constexpr double rho = 2.0;
    // const double mom = pow(-pow(lambda, rho)*log(rans[0]), 1/rho);

    const double cosT = dCos*rans[1] - 1;
    const double sinT = sqrt(1 - cosT*cosT);
    const double phi = dPhi*rans[2];
    ThreeVector pmom = {mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT};

    const double det = pow(point[1].E(), 2) + mom*mom + 2*pmom*point[1].Vec3() + pow(Constant::mN, 2);
    const double emax = Constant::mN + point[1].E() - sqrt(det);
    // const double energy = -emax/5*log(rans[3]);
    const double energy = emax*rans[3] - 1e-8;

    point[HadronIdx()] = {Constant::mN - energy, mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT};
    FourVector tmp1 = {point[HadronIdx()].Vec3(), Constant::mN - emax};
    FourVector tmp2 = {point[HadronIdx()].Vec3(), Constant::mN};
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
    spdlog::trace("  mom = {}", mom);
    spdlog::trace("  emax = {}", emax);
    spdlog::trace("  emax2 = {}", - point[1].E() - sqrt(det));
    spdlog::trace("  s = {}", (point[0] + point[1]).M2());
    spdlog::trace("  s_min = {}", (tmp1 + point[1]).M2());
    spdlog::trace("  s_max = {}", (tmp2 + point[1]).M2());
    spdlog::trace("  energy = {}", energy);
}

double QESpectralMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    constexpr double dp = 800;
    rans[0] = point[HadronIdx()].P()/dp;
    // constexpr double lambda = 200;
    // constexpr double rho = 2.0;
    // rans[0] = exp(-pow(point[HadronIdx()].P()/lambda, rho));
    // const double dp = 1.0/(rho*rans[0]*pow(point[HadronIdx()].P()/lambda, rho)/point[HadronIdx()].P());


    rans[1] = (point[HadronIdx()].CosTheta()+1)/dCos;
    rans[2] = point[HadronIdx()].Phi()/dPhi;

    const double det = pow(point[1].E(), 2) + point[0].P2() 
                     + 2*point[0].Vec3()*point[1].Vec3() + pow(Constant::mN, 2);
    const double emax = Constant::mN + point[1].E() - sqrt(det);
    const double energy = Constant::mN - point[HadronIdx()].E();
    // rans[3] = exp(-energy/(emax/5));
    // const double dE = (emax/5)/rans[3];
    rans[3] = (energy + 1e-8)/emax;
    const double dE = emax;
    double wgt = 1.0/point[0].P2()/dp/dCos/dPhi/dE;
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
    spdlog::trace("  Weight: {}", wgt);

    return wgt;
}
