#include "nuchic/Metropolis.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Utilities.hh"

// Suppress global variable warning since this is required to register to the factory
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
REGISTER_INTERACTION(nuchic::Metropolis);
using Particles = std::vector<nuchic::Particle>;

std::vector<double> nuchic::NucleonData::XSec(ChargeMode charge,
                                              const Particle &part1,
                                              const Particle&) const {
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
        // Silence magic number warning, since this is a parameterized function
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
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

Particles nuchic::NucleonData::SinglePion(ChargeMode, const std::vector<double> &rans,
                                          const Particle &part1,
                                          const Particle &part2) {
    // Get pids for all three particles
    auto PIDs = PionCharge(rans.back(), part1.ID(), part2.ID());

    // Setup outgoing particles
    Particles outgoing = { part1, part2, Particle()};
    
    // Generate final state using Rambo
    std::vector<FourVector> momenta = {part1.Momentum(), part2.Momentum(),
                                       FourVector(), FourVector(), FourVector()};
    auto masses = {outgoing[0].Mass(), outgoing[1].Mass(), outgoing[2].Mass()};
    rambo.Nout(3);
    rambo.GeneratePoint(momenta, masses, rans);

    // Set momentum and PIDs
    for(size_t i = 0; i < outgoing.size(); ++i) {
        outgoing[i].SetPID(PIDs[i]);
        outgoing[i].SetMomentum(momenta[i+2]);
    }

    return outgoing;
}

Particles nuchic::NucleonData::DoublePion(ChargeMode, const std::vector<double> &rans,
                                          const Particle &part1,
                                          const Particle &part2) {
    // Get pids for all four particles
    auto PIDs = DoublePionCharge(rans.back(), part1.ID(), part2.ID());

    // Setup outgoing particles
    Particles outgoing = { part1, part2, Particle(), Particle()};

    // Generate final state using Rambo
    std::vector<FourVector> momenta = {part1.Momentum(), part2.Momentum(),
                                       FourVector(), FourVector(), FourVector(), FourVector()};
    auto masses = {outgoing[0].Mass(), outgoing[1].Mass(), outgoing[2].Mass(), outgoing[3].Mass()};
    rambo.Nout(4);
    rambo.GeneratePoint(momenta, masses, rans);

    // Set momentum and PIDs
    for(size_t i = 0; i < outgoing.size(); ++i) {
        outgoing[i].SetPID(PIDs[i]);
        outgoing[i].SetMomentum(momenta[i+2]);
    }

    return outgoing;
}

std::vector<nuchic::PID> nuchic::NucleonData::PionCharge(const double &ran,
                                                         const PID &pid1,
                                                         const PID &pid2) const {

    auto totCharge = static_cast<int>(ParticleInfo(pid1).Charge() + ParticleInfo(pid2).Charge());

    auto frac = fpi0[totCharge % 2]; 
    if(ran < frac) {
        return {pid1, pid2, PID::pion0()};
    } else {
        if(totCharge == 1 && ran < frac + (1-frac)/2)
            return { PID::neutron(), PID::neutron(), PID::pionp() };
        else if(totCharge == 1)
            return { PID::proton(), PID::proton(), PID::pionm() };
        else {
            // The interchange of n <-> p is not needed since we generate all of the phase space
            auto pionID = totCharge == 2 ? PID::pionp() : PID::pionm();
            return { PID::neutron(), PID::proton(), pionID };
        }
    }
}

std::vector<nuchic::PID> nuchic::NucleonData::DoublePionCharge(const double &ran,
                                                               const PID &pid1,
                                                               const PID &pid2) const {

    auto totCharge = static_cast<int>(ParticleInfo(pid1).Charge() + ParticleInfo(pid2).Charge());

    if(ran < f2pi0[0]) {
        return { pid1, pid2, PID::pion0(), PID::pion0() };
    } else if(ran < f2pi0[0] + f2pi0[1]) {
        return { pid1, pid2, PID::pionp(), PID::pionm() };
    } else {
        auto sum = f2pi0[0] + f2pi0[1];
        if(totCharge == 2) {
            if(ran < f2pi0[2]+sum) 
                return { PID::proton(), PID::neutron(), PID::pionp(), PID::pion0() };
            else 
                return { PID::neutron(), PID::neutron(), PID::pionp(), PID::pionp() };
        } else if(totCharge == 0) {
            if(ran < f2pi0[2]+sum)
                return { PID::proton(), PID::neutron(), PID::pionm(), PID::pion0() };
            else
                return { PID::proton(), PID::proton(), PID::pionm(), PID::pionm() };
        } else {
            if(ran < f2pi0[3]+sum)
                return { PID::neutron(), PID::neutron(), PID::pionp(), PID::pion0() };
            else
                return { PID::proton(), PID::proton(), PID::pionm(), PID::pion0() };
        }
    }
}

std::vector<double> nuchic::PionData::XSec(ChargeMode charge,
                                           const Particle &part1,
                                           const Particle &) const {
    switch(charge) {
        case ChargeMode::ii:
            return {XSecii(part1), 0};
        case ChargeMode::ij:
            return {XSecij(part1), XSecAbs(part1)};
        case ChargeMode::zero:
            // Silence magic number warning, since this is a parameterized function
            // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
            return {0.5*XSecii(part1)+0.5*XSecij(part1), 0.5*XSecAbs(part1)};
    }
}

Particles nuchic::PionData::Elastic(ChargeMode charge,
                                       const double &kinetic,
                                       const std::array<double, 2> &rans,
                                       const Particle &part1,
                                       const Particle &part2) const {

    // Find the angular coefficient and invert the CDF numerically
    DataTable ATable;
    DataTable BTable;

    switch(charge) {
        case ChargeMode::ii:
            ATable = Aii;
            BTable = Bii;
            break;
        case ChargeMode::ij:
            ATable = Aii;
            BTable = Bii;
            break;
        case ChargeMode::zero:
            ATable = A0;
            BTable = B0;
            break;
    }

    const auto A = Lerp(kinetic, energies, ATable);
    const auto B = Lerp(kinetic, energies, BTable);
    std::function<double(double)> func = [&](double x) {
        const double norm = 5/(2*(A+5));
        // Silence magic number warning, since this is a parameterized function
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
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

Particles nuchic::PionData::SinglePion(ChargeMode mode, const std::vector<double> &rans,
                                       const Particle &part1,
                                       const Particle &part2) {
    // Get pids for all three particles
    auto PIDs = PionCharge(mode, rans.back(), part1.ID(), part2.ID());

    // Setup outgoing particles
    Particles outgoing = { part1, part2, Particle()};

    // Generate final state using Rambo
    std::vector<FourVector> momenta = {part1.Momentum(), part2.Momentum(),
                                       FourVector(), FourVector(), FourVector()};
    auto masses = {outgoing[0].Mass(), outgoing[1].Mass(), outgoing[2].Mass()};
    rambo.Nout(3);
    rambo.GeneratePoint(momenta, masses, rans);

    // Set momentum and PIDs
    for(size_t i = 0; i < outgoing.size(); ++i) {
        outgoing[i].SetPID(PIDs[i]);
        outgoing[i].SetMomentum(momenta[i+2]);
    }

    return outgoing;
}

Particles nuchic::PionData::DoublePion(ChargeMode, const std::vector<double> &rans,
                                       const Particle &part1,
                                       const Particle &part2) {
    // Get pids for all four particles
    auto PIDs = DoublePionCharge(rans.back(), part1.ID(), part2.ID());

    // Setup outgoing particles
    Particles outgoing = { part1, part2, Particle(), Particle()};

    // Generate final state using Rambo
    std::vector<FourVector> momenta = {part1.Momentum(), part2.Momentum(),
                                       FourVector(), FourVector(), FourVector(), FourVector()};
    auto masses = {outgoing[0].Mass(), outgoing[1].Mass(), outgoing[2].Mass(), outgoing[3].Mass()};
    rambo.Nout(4);
    rambo.GeneratePoint(momenta, masses, rans);

    // Set momentum and PIDs
    for(size_t i = 0; i < outgoing.size(); ++i) {
        outgoing[i].SetPID(PIDs[i]);
        outgoing[i].SetMomentum(momenta[i+2]);
    }

    return outgoing;
}

std::vector<nuchic::PID> nuchic::PionData::PionCharge(ChargeMode mode,
                                                      const double &ran,
                                                      const PID &pid1,
                                                      const PID &pid2) const {

    auto neutralProb = mode == ChargeMode::ii ? fpi0[0] : fpi0[1]; 
    if(ran < neutralProb) 
        return { pid1, pid2, PID::pion0() };
    else {
        if(ParticleInfo(pid2).Charge() == 1)
            return { pid1, PID::neutron(), PID::pionp() };
        else 
            return { pid1, PID::proton(), PID::pionm() };
    }

}

std::vector<nuchic::PID> nuchic::PionData::DoublePionCharge(const double &ran,
                                                            const PID &pid1,
                                                            const PID &pid2) const {

    if(ran < f2pi0[0]) return { pid1, pid2, PID::pion0(), PID::pion0() };
    else if(ran < f2pi0[0] + f2pi0[1]) return { pid1, pid2, PID::pionp(), PID::pionm() };
    else {
        if(ParticleInfo(pid2).Charge() == 1)
            return { pid1, PID::neutron(), PID::pion0(), PID::pionp() };
        else
            return { pid1, PID::proton(), PID::pion0(), PID::pionm() };
    }
}

bool nuchic::PionData::IsAbsorption(const std::vector<double> &xsec,
                                    const double &ran) const {
    return ran < xsec[0]/xsec[1];
}

bool nuchic::PionData::IsInelastic(ChargeMode mode,
                                   const double &kinetic,
                                   const double &ran) const {
    switch(mode) {
        case ChargeMode::ii:
            return ran < Lerp(kinetic, energies, fiiInel);
        case ChargeMode::ij:
            return ran < Lerp(kinetic, energies, fijInel);
        case ChargeMode::zero:
            return ran < Lerp(kinetic, energies, f0Inel);
    }
}

bool nuchic::PionData::IsCEX(ChargeMode mode,
                             const double &kinetic,
                             const double &ran) const {
    switch(mode) {
        case ChargeMode::ii:
            return ran < Lerp(kinetic, energies, fiiCEX);
        case ChargeMode::ij:
            return ran < Lerp(kinetic, energies, fijCEX);
        case ChargeMode::zero:
            return ran < Lerp(kinetic, energies, f0CEX);
    }
}

bool nuchic::PionData::IsSinglePion(ChargeMode,
                                    const double &kinetic,
                                    const double &ran) const {
    return ran < Lerp(kinetic, energies, fpi);
}

double nuchic::PionData::XSecii(const Particle &pion) const {
    auto kinetic = pion.Momentum().Tk();
    auto gamma = pion.E()/pion.Mass();
    if(kinetic < energies[0]) {
        // Silence magic number warning, since this is a parameterized function
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
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
        // Silence magic number warning, since this is a parameterized function
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
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
        // Silence magic number warning, since this is a parameterized function
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
        return 16.4*(0.14+eta*eta)/eta;
    } else if(kinetic > energies.back()) {
        return sigmaabs.back();
    }

    return Lerp(kinetic, energies, sigmaabs);
}

nuchic::Metropolis::Metropolis(const std::string &mode) {
    // Set the pion production mode
    ProductionMethod production; 
    if(tolower(mode) == "probe") production = ProductionMethod::Probe;
    else if(tolower(mode) == "taget") production = ProductionMethod::Target;
    else if(tolower(mode) == "random") production = ProductionMethod::Random;
    else throw std::runtime_error(fmt::format("Invalid production option: {}", mode));

    nucleon_data = std::make_shared<nuchic::NucleonData>(production); 
    pion_data = std::make_shared<nuchic::PionData>(production); 
}

double nuchic::Metropolis::CrossSection(const Particle &part1, const Particle &part2) const {
    auto xsec = CrossSections(part1, part2);

    return std::accumulate(xsec.begin(), xsec.end(), 0);
}

std::vector<double> nuchic::Metropolis::CrossSections(const Particle &part1,
                                                      const Particle &part2) const {

    auto mode = GetMode(part1, part2);
    return GetProbe(mode.first) -> XSec(mode.second, part1, part2);
}

Particles nuchic::Metropolis::GenerateFinalState(RNG &rng,
                                                 const Particle &part1,
                                                 const Particle &part2) const {
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
            outgoing = CEX(mode, rng, kinetic, part1CM, part2CM);
        } else {
            if(IsSinglePion(mode, kinetic, rng.uniform(0.0, 1.0))) {
                outgoing = SinglePion(mode, rng, part1CM, part2CM);
            }
            outgoing = DoublePion(mode, rng, part1CM, part2CM);
        }
    } else outgoing = Elastic(mode, rng, kinetic, part1CM, part2CM);

    // Return to lab frame
    for(auto part : outgoing) part.SetMomentum(part.Momentum().Boost(boostCM));

    return outgoing;
}

nuchic::InteractionType nuchic::Metropolis::GetMode(const Particle &part1, const Particle &part2) const {
    ProbeType probe = part1.Info().IsPion() ? ProbeType::pion : ProbeType::nucleon;

    if(part1.ID() == PID::pion0())
        return {probe, nuchic::ChargeMode::zero};
    else if(part1.Info().Charge() == 2*part2.Info().Charge()-1)
        return {probe, nuchic::ChargeMode::ii};
    else
        return {probe, nuchic::ChargeMode::ij};
}

std::shared_ptr<nuchic::MetropolisData> nuchic::Metropolis::GetProbe(const ProbeType &probe) const {
    return probe == ProbeType::pion ? pion_data : nucleon_data; 
}

std::vector<double> nuchic::Metropolis::XSec(InteractionType mode, const Particle &part1, const Particle &part2) const {
    return GetProbe(mode.first) -> XSec(mode.second, part1, part2);
}

// TODO: Figure out how to handle the absorption
Particles nuchic::Metropolis::Absorption(InteractionType,
                                         RNG&,
                                         const Particle &,
                                         const Particle &) const {
    // return GetProbe(mode.first) -> Absorption(mode.second, part1, part2);
    return {};
}

Particles nuchic::Metropolis::Elastic(InteractionType mode,
                                      RNG &rng,
                                      const double &kinetic,
                                      const Particle &part1,
                                      const Particle &part2) const {
    std::array<double, 2> rans{};
    rng.generate(rans, 0.0, 1.0);
    return GetProbe(mode.first) -> Elastic(mode.second, kinetic, rans, part1, part2);
}

Particles nuchic::Metropolis::CEX(InteractionType mode,
                                  RNG &rng,
                                  const double &kinetic,
                                  const Particle &part1, 
                                  const Particle &part2) const {
    std::array<double, 2> rans{};
    rng.generate(rans, 0.0, 1.0);
    return GetProbe(mode.first) -> CEX(mode.second, kinetic, rans, part1, part2);
}

Particles nuchic::Metropolis::SinglePion(InteractionType mode,
                                         RNG &rng,
                                         const Particle &part1, 
                                         const Particle &part2) const {
    // Rambo: 4*n = 12
    // Channel select: 1
    // Position generation: 1
    // Total Random numbers = 14
    constexpr size_t ndims = 14;
    std::vector<double> rans(ndims);
    rng.generate(rans, 0.0, 1.0);
    return GetProbe(mode.first) -> SinglePion(mode.second, rans, part1, part2);
}

Particles nuchic::Metropolis::DoublePion(InteractionType mode,
                                         RNG &rng,
                                         const Particle &part1, 
                                         const Particle &part2) const {
    constexpr size_t ndims = 17;
    std::vector<double> rans(ndims);
    rng.generate(rans, 0.0, 1.0);
    return GetProbe(mode.first) -> DoublePion(mode.second, rans, part1, part2);
}

bool nuchic::Metropolis::IsAbsorption(InteractionType mode, const std::vector<double> &xsec, const double &ran) const {
    return GetProbe(mode.first) -> IsAbsorption(xsec, ran);
}

bool nuchic::Metropolis::IsInelastic(InteractionType mode, const double &kinetic, const double &ran) const {
    return GetProbe(mode.first) -> IsInelastic(mode.second, kinetic, ran);
}

bool nuchic::Metropolis::IsCEX(InteractionType mode, const double &kinetic, const double &ran) const {
    return GetProbe(mode.first) -> IsInelastic(mode.second, kinetic, ran);
}

bool nuchic::Metropolis::IsSinglePion(InteractionType mode, const double &kinetic, const double &ran) const {
    return GetProbe(mode.first) -> IsInelastic(mode.second, kinetic, ran);
}
