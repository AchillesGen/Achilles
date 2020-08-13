#include "nuchic/Metropolis.hh"
#include "nuchic/Particle.hh"

bool nuchic::Metropolis::IsInelastic(const double &ran, const double &energy, bool isSame) const {
    if(energy < energies[0]) return false;

    auto idx = static_cast<size_t>(
            std::distance(energies.begin(),
                std::lower_bound(energies.begin(), energies.end(), energy)));

    auto t = energy/energies[idx];
    auto fInelMin = isSame ? fiiInel[idx] : fijInel[idx];
    auto fInelMax = isSame ? fiiInel[idx+1] : fijInel[idx+1];

    return fInelMin + t*(fInelMax-fInelMin) < ran ? true : false;
}

double nuchic::PionInteractions::xsec(const nuchic::Particle &pion,
                                      const nuchic::Particle &nucleon) const {
    auto mode = GetMode(pion, nucleon);
    auto kinetic = pion.Momentum().Tk();

    if(kinetic < energies[0]) {
        double gamma = pion.E()/pion.Mass();
        double eta = pion.Momentum().P()/pion.Mass();
        switch(mode) {
            case nuchic::PionInteractionMode::same:
                return 3.7+286*pow(gamma-1, 3);
            case nuchic::PionInteractionMode::different:
                return 6.5+23.9*(gamma-1);
            case nuchic::PionInteractionMode::zero:
                return 16.4*(0.14+eta*eta)/eta;
        }
    } else if(kinetic > energies.back()) {
        switch(mode) {
            case nuchic::PionInteractionMode::same:
                return sigmaii.back();
            case nuchic::PionInteractionMode::different:
                return sigmaij.back();
            case nuchic::PionInteractionMode::zero:
                return (sigmaij.back() + sigmaii.back())/2;
        }
    }

    auto idx = static_cast<size_t>(
            std::distance(energies.begin(),
                std::lower_bound(energies.begin(), energies.end(), kinetic)));

    auto t = kinetic/energies[idx];
    double sigmaMin{}, sigmaMax{};

    switch(mode) {
        case nuchic::PionInteractionMode::same:
            sigmaMin = sigmaii[idx];
            sigmaMax = sigmaii[idx+1];
            break;
        case nuchic::PionInteractionMode::different:
            sigmaMin = sigmaij[idx];
            sigmaMax = sigmaij[idx+1];
            break;
        case nuchic::PionInteractionMode::zero:
            sigmaMin = (sigmaii[idx] + sigmaij[idx])/2;
            sigmaMax = (sigmaii[idx+1] + sigmaij[idx+1])/2;
            break;
    }

    return sigmaMin + t*(sigmaMax - sigmaMin);
}

nuchic::PionInteractionMode nuchic::PionInteractions::GetMode(const nuchic::Particle &pion,
                                                              const nuchic::Particle &nucleon) const {
    // Get the interaction mode
    if((pion.ID() == nuchic::PID::pionp() && nucleon.ID() == nuchic::PID::proton()) 
        || (pion.ID() == nuchic::PID::pionm() && pion.ID() == nuchic::PID::neutron())) {
        return PionInteractionMode::same;
    } else if(pion.ID() == nuchic::PID::pion0()) {
        return PionInteractionMode::zero;
    }
    return PionInteractionMode::different;
}
