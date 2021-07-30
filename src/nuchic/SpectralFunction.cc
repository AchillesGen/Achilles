#include "nuchic/SpectralFunction.hh"
#include <fstream>
#include <cmath>

using nuchic::SpectralFunction;

SpectralFunction::SpectralFunction(const std::string &filename) {
    std::ifstream data(filename);
    size_t ne{}, np{};
    data >> ne >> np;
    mom.resize(np);
    energy.resize(ne);
    spectral.resize(ne*np);
    std::vector<double> dp_p(np);
    for(size_t i = 0; i < np; ++i) {
        data >> mom[i];
        for(size_t j = 0; j < ne; ++j) {
            data >> energy[j] >> spectral[i*ne+j];
        }
    }
    data.close();

    double hp = mom[1] - mom[0];
    double he = energy[1] - energy[0];
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < ne; ++j) {
            dp_p[i] += spectral[i*ne+j]*he;
        }
    }
    for(size_t i = 0; i < np; ++i) {
        norm += mom[i]*mom[i]*dp_p[i]*4*M_PI*hp;
    }

    // Find maximums for each p
    std::vector<double> maxS(np);
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < ne; ++j) {
            if(maxS[i] < spectral[i*ne+j]) maxS[i] = spectral[i*ne+j];
        }
    }

    // Setup overestimate interpolator
    overestimate = Interp1D(mom, maxS, InterpolationType::Polynomial);
    overestimate.SetPolyOrder(1);

    // Setup spectral function interpolator
    func = Interp2D(mom, energy, spectral, InterpolationType::Polynomial);
    func.SetPolyOrder(1, 1);
}
