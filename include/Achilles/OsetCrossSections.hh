#pragma once

#include "Achilles/Particle.hh"
#include <cmath>
#include <map>

namespace achilles {

class Event;

class OsetCrossSection {
  public:
    OsetCrossSection();

    double AbsCrossSection(Event &event, size_t part1, size_t part2) const;
    std::map<std::pair<PID, PID>, double> QECrossSection(Event &event, size_t part1,
                                                         size_t part2) const;

  private:
    double DensityFraction(double density) const { return density / fNormalDensity; }
    double KEFraction(double pion_KE, double pion_mass) const { return pion_KE / pion_mass; }
    double CouplingFactor(double pion_mass) const {
        return fCouplingConstant / pion_mass / pion_mass;
    }

    inline double QuadraticFunction(const double &x, const double *a) const {
        return a[0] * x * x + a[1] * x + a[2];
    }

    double SelfEnergyAbsNN(const double pion_KE, const double pion_mass,
                           const double density) const;
    double SelfEnergyAbsNNN(const double pion_KE, const double pion_mass,
                            const double density) const;
    double SelfEnergyQE(const double pion_KE, const double pion_mass, const double density) const;
    double ReducedHalfWidth(const double nucE, const double pionE, const double pion_mass,
                            const double cms_mom, const double fermimom, const double sqrts) const;
    double PXSecCommon(const double nucE, const double pionE, const double pion_mass,
                       const double cms_mom, const double fermimom, const double sqrts,
                       const double density) const;

    double avgrelvelocity(const double pion_mom, const double pionE, const double kf,
                          const double rho) const;

    const double fConstFactor = 1.0 / 12.0 / M_PI;
    const double fCouplingConstant = 0.36 * 4.0 * M_PI;
    const double fNormalDensity = 0.17; // fm-3
    const double fNormFactor = 197.327 * 197.327 * 10.0;

    // s-wave parametrization (see sec. 3.1)
    const double ImB0 = 0.035;

    //! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
    const double fCoefCQ[3] = {-5.19, 15.35, 2.06};
    //! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
    const double fCoefCA2[3] = {1.06, -6.64, 22.66};
    //! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
    const double fCoefCA3[3] = {-13.46, 46.17, -20.34};
    //! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
    const double fCoefAlpha[3] = {0.382, -1.322, 1.466};
    //! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
    const double fCoefBeta[3] = {-0.038, 0.204, 0.613};

    // s-wave parametrization (see sec. 3.1)
    const double fCoefSigma[3] = {-0.01334, 0.06889, 0.19753};
    const double fCoefB[3] = {-0.01866, 0.06602, 0.21972};
    const double fCoefD[3] = {-0.08229, 0.37062, -0.03130};
};

} // namespace achilles
