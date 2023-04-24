#ifndef SPECTRAL_FUNCTION_HH
#define SPECTRAL_FUNCTION_HH

#include <string>
#include <vector>

#include "Achilles/Interpolation.hh"

namespace achilles {

class SpectralFunction {
  public:
    SpectralFunction(const std::string &);

    // Functions
    std::vector<double> Momentum() const { return mom; }
    std::vector<double> &Momentum() { return mom; }
    std::vector<double> Energy() const { return energy; }
    std::vector<double> &Energy() { return energy; }
    double MinMomentum() const { return mom.front(); }
    double MaxMomentum() const { return mom.back(); }
    double MinEnergy() const { return energy.front(); }
    double MaxEnergy() const { return energy.back(); }
    double Normalization() const { return norm; }

    // Interpolators
    double operator()(double p) const { return overestimate(p); }
    double operator()(double p, double E) const;

  private:
    double norm{};
    std::vector<double> mom, energy, spectral;
    Interp1D overestimate;
    Interp2D func;
};

} // namespace achilles

#endif
