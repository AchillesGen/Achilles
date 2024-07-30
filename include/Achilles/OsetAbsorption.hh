#pragma once

#include "Achilles/Interpolation.hh"
#include "Achilles/Particle.hh"
#include <array>
#include <complex>
#include <iostream>
#include <map>
#include <vector>

namespace achilles {

class Event;
class FourVector;
class Potential;
class Random;

class OsetAbsCrossSection {
  public:
    OsetAbsCrossSection();

    inline double QuadraticFunction (const double &x, const double *a) const
    {
    return a[0] * x * x + a[1] * x + a[2];
    }

    double AbsorptionCrossSection(Event &event, size_t part1, size_t part2) const;

  private:
    const double fCouplingConstant = 0.36 * 4.0 * M_PI;
    const double fNormalDensity = 0.17; // fm-3
    const double fNormFactor = 197.327 * 197.327 * 10.0;

    // s-wave parametrization (see sec. 3.1)
    const double ImB0 = 0.035;

    //! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
    const double fCoefCQ[3] = { -5.19, 15.35,   2.06};
    //! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
    const double fCoefCA2[3] = {  1.06, -6.64,  22.66};
    //! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
    const double fCoefCA3[3] = {-13.46, 46.17, -20.34};
    //! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
    const double fCoefAlpha[3] = {0.382, -1.322, 1.466};
    //! delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
    const double fCoefBeta[3]  = {-0.038,  0.204, 0.613};

    double DeltaReduction (const double pion_E, const double pion_mass, const double momentum_cms, const double fermi_energy, const double sqrts) const;

};

} // // namespace achilles
