#include "Achilles/SpectralFunction.hh"
#include "spdlog/spdlog.h"
#include <cmath>
#include <fstream>

using achilles::SpectralFunction;

SpectralFunction::SpectralFunction(const std::string &filename) {
    std::ifstream data(filename);
    size_t ne{}, np{};
    data >> ne >> np;
    mom.resize(np);
    energy.resize(ne);
    spectral.resize(ne * np);
    std::vector<double> dp_p(np);
    for(size_t j = 0; j < np; ++j) {
        data >> mom[j];
        for(size_t i = 0; i < ne; ++i) { data >> energy[i] >> spectral[j * ne + i]; }
    }
    data.close();

    double hp = mom[1] - mom[0];
    double he = energy[1] - energy[0];
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < ne; ++j) { dp_p[i] += spectral[i * ne + j] * he; }
    }
    for(size_t i = 0; i < np; ++i) { norm += mom[i] * mom[i] * dp_p[i] * 4 * M_PI * hp; }
    spdlog::debug("Spectral function normalization: {}", norm);

    // Find maximums for each p
    std::vector<double> maxS(np);
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < ne; ++j) {
            if(maxS[i] < spectral[i * ne + j]) maxS[i] = spectral[i * ne + j];
        }
    }

    // Setup overestimate interpolator
    overestimate = Interp1D(mom, maxS, InterpolationType::Polynomial);
    overestimate.SetPolyOrder(1);

    // Setup spectral function interpolator
    func = Interp2D(mom, energy, spectral, InterpolationType::Polynomial);
    func.SetPolyOrder(3, 1);
}

double SpectralFunction::operator()(double p, double E) const {
    if(p < mom.front() || p > mom.back() || E < energy.front() || E > energy.back()) return 0;

    auto result = func(p, E);
    return result > 0 ? result : 0;
}
