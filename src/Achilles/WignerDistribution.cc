#include "Achilles/WignerDistribution.hh"
#include "Achilles/System.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Utilities.hh"
#include "spdlog/spdlog.h"
#include <set>
#include <fstream>

using achilles::WignerDistribution;

WignerDistribution::WignerDistribution(const std::string &filename) {
    spdlog::debug("Reading Wigner Distribution from file: {}", filename);
    std::ifstream data(Filesystem::FindFile(filename, "WignerDistribution"));
    size_t nr{}, np{};
    data >> nr >> np;
    mom.resize(np);
    radius.resize(nr);
    wigner.resize(nr * np);
    dr_r.resize(nr);    
    std::vector<double> dp_p(np);

    double k,r,w,e;
    std::set<double> radius_set;
    std::set<double> mom_set;

    for(size_t j = 0; j < nr; ++j) {
        for(size_t i = 0; i < np; ++i) {
            data >> r >> k >> w >> e;
            radius_set.insert(r);
            mom_set.insert(k * Constant::HBARC);
            wigner[j * np + i] = (w / pow(Constant::HBARC,3));
        }
    }

    data.close();

    std::copy(radius_set.begin(), radius_set.end(), radius.begin());
    std::copy(mom_set.begin(), mom_set.end(), mom.begin());

    double hp = mom[1] - mom[0];
    double hr = radius[1] - radius[0];

    // Compute n(k)
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < nr; ++j) { 
            dp_p[i] += wigner[j * np + i] * radius[j] * radius[j] * 4 * M_PI * hr; 
        }
    }

    // Compute rho(r)
    for(size_t j = 0; j < nr; ++j) {
        for(size_t i = 0; i < np; ++i) { dr_r[j] += std::abs(wigner[j * np + i]) * mom[i] * mom[i] * 4 * M_PI * hp / pow(2 * M_PI,3);}
    }

    for(size_t i = 0; i < np; ++i) { norm += mom[i] * mom[i] * dp_p[i] * 4 * M_PI * hp / pow(2 * M_PI,3); }


    // Setup momentum distribution interpolator
    rho_func = Interp1D(radius, dr_r, InterpolationType::Polynomial);
    rho_func.SetPolyOrder(3);

    // Setup spectral function interpolator
    func = Interp2D(mom, radius, wigner, InterpolationType::Polynomial);
    func.SetPolyOrder(3, 1);
}

double WignerDistribution::operator()(double p, double r) const {
    if(p < mom.front() || p > mom.back() || r < radius.front() || r > radius.back()) return 0;

    auto result = func(p, r) / norm;
    return result;
    //return result > 0 ? result : 0;
}

double WignerDistribution::operator()(double r) const {
    if(r < radius.front() || r > radius.back()) return 0;

    auto result = rho_func(r) / norm;
    return result > 0 ? result : 0;
}

double WignerDistribution::MaxWeight(double r) const {
    auto absfunc = [&](double p) {return -std::abs(this->operator()(p,r));};
    Brent brent(absfunc);
    return -absfunc(brent.Minimize(mom.front(), mom.back()));
}
