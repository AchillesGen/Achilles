#include "Achilles/SpectralFunction.hh"
#include "spdlog/spdlog.h"
#include <fstream>
#include <cmath>

using achilles::SpectralFunction;

SpectralFunction::SpectralFunction(const std::string &filename) {
    spdlog::debug("Reading spectral function from: {}", filename);
    std::ifstream data(filename);
    size_t ne{}, np{};
    
    data >> ne >> np;
    mom.resize(np);
    energy.resize(ne);
    spectral.resize(ne*np);
    spectral_mf.resize(ne*np);
    spectral_corr.resize(ne*np);
    dp_p_mf.resize(np);
    dp_p_corr.resize(np);
    dp_p.resize(np);

    // Quick and dirty cut off because 
    // No QMC SF provided
    double cutoff = 300.; //MF-Corr cutoff

    for(size_t j = 0; j < np; ++j) {
        data >> mom[j];
        for(size_t i = 0; i < ne; ++i) {
            data >> energy[i] >> spectral[j*ne+i];
            //Fill MF and Corr SF's
            mom[j] < cutoff ? spectral_mf[j*ne +i] = spectral[j*ne+i] : spectral_corr[j*ne +i] = spectral[j*ne+i];

        }
    }
    data.close();
    
    // TODO
    // members to add
    // meanfield + correlation SF

    
    
    /*
    data >> ne;
    data >> np;
    ne=400;
    np=200;	    
    mom.resize(np);
    energy.resize(ne);
    spectral.resize(ne*np);
    std::vector<double> dp_p(np);

    for (size_t j = 0; j < np; ++j) {
        for (size_t i = 0; i < ne; ++i) {
            data >> mom[j] >> energy[i] >> spectral[j*ne+i];
        }
    }


    data.close();
    */

    double hp = mom[1] - mom[0];
    double he = energy[1] - energy[0];
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < ne; ++j) {
            dp_p[i] += spectral[i*ne+j]*he;
            dp_p_mf[i] += spectral_mf[i*ne+j]*he;
            dp_p_corr[i] += spectral_corr[i*ne+j]*he;
        }
    }
    for(size_t i = 0; i < np; ++i) {
        norm += mom[i]*mom[i]*dp_p[i]*4*M_PI*hp;
    }
    spdlog::debug("Spectral function normalization: {}", norm);

    // Find maximums for each p
    std::vector<double> maxS(np);
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < ne; ++j) {
            if(maxS[i] < spectral[i*ne+j]) maxS[i] = spectral[i*ne+j];
        }
    }

    // Setup momentum distribution interpolator
    momentum_distribution = Interp1D(mom, dp_p, InterpolationType::Polynomial);
    momentum_distribution.SetPolyOrder(3);

    // Setup spectral function interpolator
    func = Interp2D(mom, energy, spectral, InterpolationType::Polynomial);
    func.SetPolyOrder(3, 1);
}

double SpectralFunction::operator()(double p, double E) const {
    if(p < mom.front() || p > mom.back() || E < energy.front() || E > energy.back())
        return 0;

    auto result = func(p, E);
    return result > 0 ? result : 0;
}

double SpectralFunction::operator()(double p) const {
    if(p < mom.front() || p > mom.back())
        return 0;


    auto result = momentum_distribution(p);
    return result > 0 ? result : 0;

}
