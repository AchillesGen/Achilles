#ifndef WIGNER_DISTRIBUTION_HH
#define WIGNER_DISTRIBUTION_HH

#include <string>
#include <vector>

#include "Achilles/Interpolation.hh"

namespace achilles {

class WignerDistribution {
  public:
    WignerDistribution(const std::string &);

    // Functions
    std::vector<double> Momentum() const { return mom; }
    std::vector<double> &Momentum() { return mom; }
    std::vector<double> Radius() const { return radius; }
    std::vector<double> &Radius() { return radius; }
    double MinMomentum() const { return mom.front(); }
    double MaxMomentum() const { return mom.back(); }
    double MinRadius() const { return radius.front(); }
    double MaxRadius() const { return radius.back(); }
    double Normalization() const { return norm; }

    // Interpolators
    double operator()(double p) const;
    double operator()(double p, double r) const;

  private:
    double norm{};
    std::vector<double> mom, radius, wigner, dp_p;
    Interp1D mom_func;
    Interp2D func;
};

} // namespace achilles

#endif
