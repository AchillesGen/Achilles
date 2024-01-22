#ifndef SPECTRAL_FUNCTION_HH
#define SPECTRAL_FUNCTION_HH

#include <string>
#include <vector>

#include "Achilles/Interpolation.hh"

namespace achilles {

class SpectralFunction {
    public:
        SpectralFunction(const std::string&);

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
        double operator()(double p) const;
        double operator()(double p, double E) const;

    private:
        double norm{};
        std::vector<double> mom, energy, spectral, dp_p;
        std::vector<double> spectral_mf, spectral_corr;
        std::vector<double> dp_p_mf, dp_p_corr;
        Interp1D momentum_distribution;
        Interp2D func;
};

}

#endif
