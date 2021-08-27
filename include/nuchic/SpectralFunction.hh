#ifndef SPECTRAL_FUNCTION_HH
#define SPECTRAL_FUNCTION_HH

#include "nuchic/Interpolation.hh"
#include <string>

namespace nuchic {

class SpectralFunction {
    public:
        SpectralFunction(const std::string&);

        double operator()(double) const;

    private:
        Interp1D m_spectral;
};

}

#endif
