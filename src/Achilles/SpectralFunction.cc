#include "Achilles/SpectralFunction.hh"
#include "spdlog/spdlog.h"
#include <fstream>

using achilles::SpectralFunction;

SpectralFunction::SpectralFunction(const std::string &filename) {
    std::ifstream data(filename);
    size_t ne{}, np{};
    data >> ne >> np;
    mom.resize(np);
    energy.resize(ne);
    spectral.resize(ne * np);
    dp_p.resize(np);
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

    // Setup momentum distribution interpolator
    mom_func = Interp1D(mom, dp_p, InterpolationType::Polynomial);
    mom_func.SetPolyOrder(3);

    // Setup spectral function interpolator
    func = Interp2D(mom, energy, spectral, InterpolationType::Polynomial);
    func.SetPolyOrder(3, 1);
}

double SpectralFunction::operator()(double p, double E) const {
    if(p < mom.front() || p > mom.back() || E < energy.front() || E > energy.back()) return 0;

    auto result = func(p, E) / norm;
    return result > 0 ? result : 0;
}

double SpectralFunction::operator()(double p) const {
    if(p < mom.front() || p > mom.back()) return 0;

    auto result = mom_func(p) / norm;
    return result > 0 ? result : 0;
}
