#include "nuchic/HadronicMapper.hh"
#include "nuchic/Constants.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"
#include "spdlog/spdlog.h"

using nuchic::QESpectralMapper;

void QESpectralMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    // Generate inital nucleon state
    // const double mom = dp*rans[0];
    constexpr double lambda = 200;
    constexpr double rho = 1.5;
    const double mom = pow(-pow(lambda, rho)*log(rans[0]), 1/rho);
    const double cosT = dCos*rans[1] - 1;
    const double sinT = sqrt(1 - cosT*cosT);
    const double phi = dPhi*rans[2];
    ThreeVector pmom = {mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT};

    const double det = pow(point[1].E(), 2) + mom*mom + 2*pmom*point[1].Vec3();
    const double emax = Constant::mN + point[1].E() - sqrt(det);
    const double energy = -emax/10*log(rans[3]);
    // const double energy = emax*rans[3] - 1e-8;
    // std::cout << "QEMapper: " << emax << " " << Constant::mN + point[1].E() << " " << sqrt(det) << " " << rans[3] << " ";

    point[HadronIdx()] = {Constant::mN - energy, mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT};
    FourVector tmp1 = {point[HadronIdx()].Vec3(), Constant::mN - emax};
    FourVector tmp2 = {point[HadronIdx()].Vec3(), Constant::mN};
    // std::cout << (point[0] + point[1]).M2() << " " << (tmp1 + point[1]).M2() << std::endl;
    spdlog::trace("QESpectralMapper::GeneratePoint:");
    spdlog::trace("  rans = [{}]", fmt::join(rans.begin(), rans.end(), ", "));
    spdlog::trace("  mom = {}", mom);
    spdlog::trace("  emax = {}", emax);
    spdlog::trace("  emax2 = {}", - point[1].E() - sqrt(det));
    spdlog::trace("  s = {}", (point[0] + point[1]).M2());
    spdlog::trace("  s_max = {}", (tmp1 + point[1]).M2());
    spdlog::trace("  s_min = {}", (tmp2 + point[1]).M2());
    spdlog::trace("  energy = {}", energy);
}

double QESpectralMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    constexpr double lambda = 300;
    constexpr double rho = 1.5;
    rans[0] = exp(-pow(point[HadronIdx()].P()/lambda, rho));
    const double dp = rho*rans[0]*pow(point[HadronIdx()].P()/lambda, rho)/point[HadronIdx()].P();

    rans[1] = (point[HadronIdx()].CosTheta()+1)/dCos;
    rans[2] = point[HadronIdx()].Phi()/dPhi;

    const double det = pow(point[1].E(), 2) + point[0].P2() + 2*point[0].Vec3()*point[1].Vec3();
    const double emax = Constant::mN + point[1].E() - sqrt(det);
    const double energy = Constant::mN - point[HadronIdx()].E();
    rans[3] = exp(-energy/(emax/10));
    const double dE = rans[3]/(emax/10);

    return 1.0/point[HadronIdx()].P2()*dp/dCos/dPhi*dE;
}
