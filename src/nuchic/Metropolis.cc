#include "nuchic/Metropolis.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Utilities.hh"

REGISTER_INTERACTION(nuchic::Metropolis);
using Particles = std::vector<nuchic::Particle>;

std::vector<double> nuchic::NucleonData::XSec(ChargeMode charge,
                                              const Particle &part1,
                                              const Particle &part2) const {
    auto sigma = (charge == ChargeMode::ii) ? sigmaii : sigmaij;
    auto Tk = part1.Momentum().Tk();
    return { Lerp(Tk, energies, sigma) };
}

bool nuchic::NucleonData::IsInelastic(ChargeMode charge,
                                      const double &kinetic,
                                      const double &ran) const {
    auto fInel = (charge == ChargeMode::ii) ? fiiInel : fijInel;
    return Lerp(kinetic, energies, fInel) < ran;
}

bool nuchic::NucleonData::IsSinglePion(ChargeMode,
                                       const double &kinetic,
                                       const double &ran) const {
    return Lerp(kinetic, energies, fpi) < ran; 
}

Particles nuchic::NucleonData::Elastic(ChargeMode charge,
                                       const double &kinetic,
                                       const std::array<double, 2> &rans,
                                       const Particle &part1,
                                       const Particle &part2) const {

    // Find the angular coefficient and invert the CDF numerically
    const auto ATable = (charge == ChargeMode::ii) ? Aii : Aij;
    const auto BTable = (charge == ChargeMode::ii) ? Bii : Bij;
    const auto A = Lerp(kinetic, energies, ATable);
    const auto B = Lerp(kinetic, energies, BTable);
    std::function<double(double)> func = [&](double x) {
        const double norm = 5/(2*(A+5));
        return norm*(A/5*pow(x, 5)+B/4*pow(x, 4)+x+A/5-B/4+1)-rans[0];
    };
    auto brent = Brent(func);

    double pcm = part1.Momentum().P();
    double costheta = brent.CalcRoot(-1, 1);
    double sintheta = std::sqrt(1-costheta*costheta);
    double phi = 2*M_PI*rans[1];

    // Setup the outgoing in the CM frame
    Particles outgoing = { part1, part2 };
    outgoing[0].SetMomentum({pcm*sintheta*cos(phi),
                             pcm*sintheta*sin(phi),
                             pcm*costheta, part1.E()});
    outgoing[1].SetMomentum({-pcm*sintheta*cos(phi),
                             -pcm*sintheta*sin(phi),
                             -pcm*costheta, part1.E()});

    return outgoing;
}

Particles nuchic::NucleonData::SinglePion(ChargeMode charge, const std::array<double, 5> &rans,
                                          const Particle &part1,
                                          const Particle &part2) const {
    // Evenly split momentum over all three particles
    double pcm = part1.Momentum().P()*2/3; 
    double ecm = part1.E()*2/3;

    // Angles for first particle
    double costheta1 = 2*rans[0] - 1;
    double sintheta1 = std::sqrt(1-costheta1*costheta1);
    double phi1 = 2*M_PI*rans[1];

    // Angles for second particle
    double costheta2 = 2*rans[2] - 1;
    double sintheta2 = std::sqrt(1-costheta2*costheta2);
    double phi2 = 2*M_PI*rans[3];

    // Get pids for all three particles
    int totCharge = static_cast<int>(part1.Info().Charge() + part2.Info().Charge());
    int pionCharge = PionCharge(rans[4], totCharge);
}

int nuchic::NucleonData::PionCharge(const double &ran, int totCharge) const {
    auto frac = fpi0[totCharge % 2]; 
    if(ran < frac) {
        return 0;
    } else {
        if(totCharge == 1 && ran < frac + (1-frac)/2)
            return 1;
        else if(totCharge == 1)
            return -1;
        else
            return totCharge-1;
    }
}

std::pair<int, int> nuchic::NucleonData::DoublePionCharge(const double &ran,
                                                          int totCharge) const {
    if(ran < f2pi0[0]) {
        return { 0, 0 };
    } else if(ran < f2pi0[0] + f2pi0[1]) {
        return { 1, -1 };
    } else {
        ran -= f2pi0[0] + f2pi0[1];
        if(totCharge == 2) {
            if(ran < f2pi0[2]) return { 1, 0 };
            else return { 1, 1 };
        } else if(totCharge == 0) {
            if(ran < f2pi0[2]) return { -1, 0 };
            else return { -1, -1 };
        } else {
            if(ran < f2pi0[3]) return { 1, 0 };
            else return { -1 , 0 };
        }
    }
}

std::vector<double> nuchic::PionData::XSec(int mode, const Particle &pion) const {
    switch(mode) {
        case 2:
            return {XSecii(pion), 0};
        case 1:
            return {XSecij(pion), XSecAbs(pion)};
        case 0:
            return {0.5*XSecii(pion)+0.5*XSecij(pion),
                    0.5*XSecAbs(pion)};
    }
}

bool nuchic::PionData::IsAbsorption(const double &ran,
                                    const std::vector<double> &xsec) const {
    return ran < xsec[0]/xsec[1];
}

bool nuchic::PionData::IsInelastic(const double &ran,
                                   const nuchic::Particle &pion,
                                   const nuchic::Particle &nucleon) const {
    auto mode = GetMode(pion, nucleon);
    auto kinetic = pion.Momentum().Tk();
    switch(mode) {
        case nuchic::PionInteractionMode::same:
            return ran < Lerp(kinetic, energies, fiiInel);
        case nuchic::PionInteractionMode::different:
            return ran < Lerp(kinetic, energies, fijInel);
        case nuchic::PionInteractionMode::zero:
            return ran < Lerp(kinetic, energies, f0Inel);
    }
}

bool nuchic::PionData::IsCEX(const double &ran,
                             const nuchic::Particle &pion,
                             const nuchic::Particle &nucleon) const {
    auto mode = GetMode(pion, nucleon);
    auto kinetic = pion.Momentum().Tk();
    switch(mode) {
        case nuchic::PionInteractionMode::same:
            return ran < Lerp(kinetic, energies, fiiCEX);
        case nuchic::PionInteractionMode::different:
            return ran < Lerp(kinetic, energies, fijCEX);
        case nuchic::PionInteractionMode::zero:
            return ran < Lerp(kinetic, energies, f0CEX);
    }
}

bool nuchic::PionData::IsSinglePion(const double &ran,
                                    const nuchic::Particle &pion,
                                    const nuchic::Particle&) const {
    auto kinetic = pion.Momentum().Tk();
    return ran < Lerp(kinetic, energies, fpi);
}

double nuchic::PionData::XSecii(const Particle &pion) const {
    auto kinetic = pion.Momentum().Tk();
    auto gamma = pion.E()/pion.Mass();
    if(kinetic < energies[0]) {
        return 3.7+286*pow(gamma-1, 3);
    } else if(kinetic > energies.back()) {
        return sigmaii.back();
    }

    return Lerp(kinetic, energies, sigmaii);
}

double nuchic::PionData::XSecij(const Particle &pion) const {
    auto kinetic = pion.Momentum().Tk();
    auto gamma = pion.E()/pion.Mass();
    if(kinetic < energies[0]) {
        return 6.5+23.9*(gamma-1);
    } else if(kinetic > energies.back()) {
        return sigmaij.back();
    }

    return Lerp(kinetic, energies, sigmaij);
}

double nuchic::PionData::XSecAbs(const Particle &pion) const {
    auto kinetic = pion.Momentum().Tk();
    double eta = pion.Momentum().P()/pion.Mass();
    if(kinetic < energies[0]) {
        return 16.4*(0.14+eta*eta)/eta;
    } else if(kinetic > energies.back()) {
        return sigmaabs.back();
    }

    return Lerp(kinetic, energies, sigmaabs);
}

nuchic::Metropolis::Metropolis(const std::string&) {
    nucleon_data = std::make_shared<nuchic::NucleonData>(); 
    pion_data = std::make_shared<nuchic::PionData>(); 
}

double nuchic::Metropolis::CrossSection(const Particle &part1, const Particle &part2) const {
    auto xsec = CrossSections(part1, part2);

    return std::accumulate(xsec.begin(), xsec.end(), 0);
}

std::vector<double> nuchic::Metropolis::CrossSections(const Particle &part1,
                                                      const Particle &part2) const {
    if(part1.ID() == PID::neutron() || part1.ID() == PID::proton()) {
        return nucleon_data -> XSec(part1, part2);
    } else {
        return pion_data -> XSec(part1, part2);
    }
}

Particles nuchic::Metropolis::GenerateFinalState(randutils::mt19937_rng &rng,
                                                 const Particle &part1,
                                                 const Particle &part2) {
    auto mode = GetMode(part1, part2);
    auto kinetic = part1.Momentum().Tk();

    // Boost to the CM frame
    auto part1CM = part1;
    auto part2CM = part2;
    auto boostCM = (part1.Momentum() + part2.Momentum()).BoostVector();
    part1CM.SetMomentum(part1CM.Momentum().Boost(-boostCM));
    part2CM.SetMomentum(part2CM.Momentum().Boost(-boostCM));

    Particles outgoing{};

    // TODO: Figure out how to handle absorption
    // if(IsAbsorption(mode, kinetic, rng.uniform(0.0, 1.0))) 
    //     outgoing =  Absorption(mode, part1CM, part2CM);
    if(IsInelastic(mode, kinetic, rng.uniform(0.0, 1.0))) {
        if(IsCEX(mode, kinetic, rng.uniform(0.0, 1.0))) {
            outgoing = CEX(mode, part1CM, part2CM);
        } else {
            if(IsSinglePion(mode, kinetic, rng.uniform(0.0, 1.0))) {
                outgoing = SinglePion(mode, part1CM, part2CM);
            }
            outgoing = DoublePion(mode, part1CM, part2CM);
        }
    } else outgoing = Elastic(mode, kinetic, part1CM, part2CM);

    // Return to lab frame
    for(auto part : outgoing) part.SetMomentum(part.Momentum().Boost(boostCM));

    return outgoing;
}

nuchic::ChargeMode nuchic::Metropolis::GetMode(const Particle &part1, const Particle &part2) {
    auto pid1 = part1.ID(), pid2 = part2.ID();

    if(abs(pid1) < 1000) current_data = pion_data;
    else current_data = nucleon_data;

    if(pid1 == PID::pion0())
        return nuchic::ChargeMode::zero;
    else if((pid1 == PID::pionp() && pid2 == PID::proton())
            || (pid1 == PID::pionm() && pid2 == PID::neutron())
            || (pid1 == pid2))
        return nuchic::ChargeMode::ii;
    else
        return nuchic::ChargeMode::ij;
}

std::vector<double> nuchic::Metropolis::XSec(ChargeMode mode, const Particle &part1) const {
    return current_data -> XSec(mode, part1);
}

// TODO: Figure out how to handle the absorption
Particles nuchic::Metropolis::Absorption(ChargeMode mode,
                                         const Particle &part1,
                                         const Particle &part2) const {
}

Particles nuchic::Metropolis::Elastic(ChargeMode mode,
                                      const double &kinetic,
                                      const Particle &part1,
                                      const Particle &part2) const {
    return current_data -> Elastic(mode, kinetic, part1, part2);
}

Particles nuchic::Metropolis::CEX(ChargeMode mode,
                                  const Particle &part1, 
                                  const Particle &part2) const {
    return current_data -> CEX(mode, part1, part2);
}

Particles nuchic::Metropolis::SinglePion(ChargeMode mode,
                                         const Particle &part1, 
                                         const Particle &part2) const {
    return current_data -> SinglePion(mode, part1, part2);
}

Particles nuchic::Metropolis::DoublePion(ChargeMode mode,
                                         const Particle &part1, 
                                         const Particle &part2) const {
    return current_data -> DoublePion(mode, part1, part2);
}

bool nuchic::Metropolis::IsAbsorption(ChargeMode mode, const double &kinetic, const double &ran) const {
    return current_data -> IsAbsorption(mode, kinetic, ran);
}

bool nuchic::Metropolis::IsInelastic(ChargeMode mode, const double &kinetic, const double &ran) const {
    return current_data -> IsInelastic(mode, kinetic, ran);
}

bool nuchic::Metropolis::IsCEX(ChargeMode mode, const double &kinetic, const double &ran) const {
    return current_data -> IsInelastic(mode, kinetic, ran);
}

bool nuchic::Metropolis::IsSinglePion(ChargeMode mode, const double &kinetic, const double &ran) const {
    return current_data -> IsInelastic(mode, kinetic, ran);
}
