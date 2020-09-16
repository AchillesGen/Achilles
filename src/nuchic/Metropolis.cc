#include "nuchic/Metropolis.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Utilities.hh"

using DataTable = std::array<double, 8>;

// This is required to make sure the symbols are defined

// NucleonData static members
constexpr DataTable nuchic::NucleonData::energies;
constexpr DataTable nuchic::NucleonData::sigmaii;
constexpr DataTable nuchic::NucleonData::sigmaij;
constexpr DataTable nuchic::NucleonData::fiiInel;
constexpr DataTable nuchic::NucleonData::fijInel;
constexpr DataTable nuchic::NucleonData::fpi;
constexpr DataTable nuchic::NucleonData::Aii;
constexpr DataTable nuchic::NucleonData::Bii;
constexpr DataTable nuchic::NucleonData::Aij;
constexpr DataTable nuchic::NucleonData::Bij;
constexpr std::array<double, 2> nuchic::NucleonData::fpi0;
constexpr std::array<double, 4> nuchic::NucleonData::f2pi0;

// PionData static members
constexpr DataTable nuchic::PionData::energies;
constexpr DataTable nuchic::PionData::sigmaii;
constexpr DataTable nuchic::PionData::sigmaij;
constexpr DataTable nuchic::PionData::sigmaabs;
constexpr DataTable nuchic::PionData::fiiInel;
constexpr DataTable nuchic::PionData::fiiCEX;
constexpr DataTable nuchic::PionData::fijInel;
constexpr DataTable nuchic::PionData::fijCEX;
constexpr DataTable nuchic::PionData::f0Inel;
constexpr DataTable nuchic::PionData::f0CEX;
constexpr DataTable nuchic::PionData::fpi;
constexpr DataTable nuchic::PionData::Aii;
constexpr DataTable nuchic::PionData::Bii;
constexpr DataTable nuchic::PionData::Aij;
constexpr DataTable nuchic::PionData::Bij;
constexpr DataTable nuchic::PionData::A0;
constexpr DataTable nuchic::PionData::B0;
constexpr std::array<double, 2> nuchic::PionData::fpi0;
constexpr std::array<double, 2> nuchic::PionData::f2pi0;

// Suppress global variable warning since this is required to register to the factory
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
REGISTER_INTERACTION(nuchic::Metropolis);
using Particles = std::vector<nuchic::Particle>;

nuchic::ThreeVector nuchic::MetropolisData::GetPosition(ProductionMethod mode,
                                                        const double &ran,
                                                        const ThreeVector &r1,
                                                        const ThreeVector &r2) const {
    ThreeVector position{};
    switch(mode) {
        case ProductionMethod::Probe:
            position = r1;
            break;
        case ProductionMethod::Target:
            position = r2;
            break;
        case ProductionMethod::Random:
            position = (r2-r1)*ran + r1;
            break;
    }

    return position;
}

std::vector<double> nuchic::NucleonData::XSec(ChargeMode charge,
                                              const Particle &part1,
                                              const Particle &part2) const {
    auto sigma = (charge == ChargeMode::ii) ? sigmaii : sigmaij;
    FourVector momentum = part1.Momentum().Boost(-part2.Momentum().BoostVector());
    auto Tk = momentum.Tk();

    if(Tk < energies[0]) {
        auto beta = momentum.BoostVector().P();
        if(charge == ChargeMode::ii) {
            // Silence magic number warning, since this is a parameterized function
            // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
            return {10.63/beta/beta - 29.92/beta + 42.9};
        } else {
            // Silence magic number warning, since this is a parameterized function
            // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
            return {34.10/beta/beta - 82.20/beta + 82.2};
        }
    }
    return { Lerp(Tk, energies, sigma) };
}

bool nuchic::NucleonData::IsInelastic(ChargeMode charge,
                                      const double &kinetic,
                                      const double &ran) const {
    auto fInel = (charge == ChargeMode::ii) ? fiiInel : fijInel;
    return Lerp(kinetic, energies, fInel) > ran;
}

bool nuchic::NucleonData::IsSinglePion(ChargeMode,
                                       const double &kinetic,
                                       const double &ran) const {
    return Lerp(kinetic, energies, fpi) > ran; 
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
    double costheta = 0;
    try {
        costheta = brent.CalcRoot(-1, 1);
    } catch (const std::domain_error &e) {
        costheta = pow(2*rans[0]-1, 1.0/5);
    }
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
    auto masses = {part1.Mass(), part2.Mass(),
                   outgoing[0].Mass(), outgoing[1].Mass(),
                   ParticleInfo(PID::pionp()).Mass()};
    rambo.Nout(3);
    rambo.GeneratePoint(momenta, masses, rans);

    // Set position for pion
    outgoing[2].SetPosition(GetPosition(production, rans[rans.size()-2],
                                        part1.Position(), part2.Position()));

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
    auto masses = {part1.Mass(), part2.Mass(),
                   outgoing[0].Mass(), outgoing[1].Mass(),
                   ParticleInfo(PID::pionp()).Mass(), ParticleInfo(PID::pionp()).Mass()};
    rambo.Nout(4);
    rambo.GeneratePoint(momenta, masses, rans);

    // Set position for pion
    outgoing[2].SetPosition(GetPosition(production, rans[rans.size()-2],
                                        part1.Position(), part2.Position()));
    outgoing[3].SetPosition(GetPosition(production, rans[rans.size()-3],
                                        part1.Position(), part2.Position()));

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

    auto totCharge = static_cast<size_t>(ParticleInfo(pid1).Charge() + ParticleInfo(pid2).Charge());

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

    auto totCharge = static_cast<size_t>(ParticleInfo(pid1).Charge() + ParticleInfo(pid2).Charge());

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
                                           const Particle &part2) const {

    Particle restFrame = part1;
    restFrame.Momentum().Boost(-part2.Momentum().BoostVector());

    switch(charge) {
        case ChargeMode::ii:
            return { XSecii(restFrame), 0 };
        case ChargeMode::ij:
            return { XSecij(restFrame), XSecAbs(restFrame) };
        case ChargeMode::zero:
            // Silence magic number warning, since this is a parameterized function
            // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
            return { 0.5*XSecii(restFrame)+0.5*XSecij(restFrame), 0.5*XSecAbs(restFrame) };
    }

    return {};
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
    double costheta = 0;
    try {
        costheta = brent.CalcRoot(-1, 1);
    } catch (const std::domain_error &e) {
        costheta = pow(2*rans[0]-1, 1.0/5);
    }
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

Particles nuchic::PionData::CEX(ChargeMode charge,
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
    double costheta = 0;
    try {
        costheta = brent.CalcRoot(-1, 1);
    } catch (const std::domain_error &e) {
        costheta = pow(2*rans[0]-1, 1.0/5);
    }
    double sintheta = std::sqrt(1-costheta*costheta);
    double phi = 2*M_PI*rans[1];

    // Setup the outgoing in the CM frame
    Particles outgoing = { part1, part2 };
    if(part2.ID() == PID::proton()) {
        outgoing[0].SetPID(PID::pionp());
        outgoing[1].SetPID(PID::neutron());
    } else {
        outgoing[0].SetPID(PID::pionm());
        outgoing[1].SetPID(PID::proton());
    }
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
    auto masses = {part1.Mass(), part2.Mass(),
                   outgoing[0].Mass(), outgoing[1].Mass(),
                   ParticleInfo(PID::pionp()).Mass()};
    rambo.Nout(3);
    rambo.GeneratePoint(momenta, masses, rans);

    // Set position for pion
    outgoing[2].SetPosition(GetPosition(production, rans[rans.size()-2],
                                        part1.Position(), part2.Position()));

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
    auto masses = {part1.Mass(), part2.Mass(),
                   outgoing[0].Mass(), outgoing[1].Mass(),
                   ParticleInfo(PID::pionp()).Mass(), ParticleInfo(PID::pionp()).Mass()};
    rambo.Nout(4);
    rambo.GeneratePoint(momenta, masses, rans);

    // Set position for pion
    outgoing[2].SetPosition(GetPosition(production, rans[rans.size()-2],
                                        part1.Position(), part2.Position()));
    outgoing[3].SetPosition(GetPosition(production, rans[rans.size()-3],
                                        part1.Position(), part2.Position()));

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

    return false;
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

    return false;
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

    auto energy = part1CM.E() + part2CM.E() - part1.Mass() - part2.Mass();

    Particles outgoing{};

    // TODO: Figure out how to handle absorption and clean up all the if statements
    // if(IsAbsorption(mode, kinetic, rng.uniform(0.0, 1.0))) 
    //     outgoing =  Absorption(mode, part1CM, part2CM);
    if(IsInelastic(mode, kinetic, rng.uniform(0.0, 1.0))) {
        if(IsCEX(mode, kinetic, rng.uniform(0.0, 1.0))) {
            outgoing = CEX(mode, rng, kinetic, part1CM, part2CM);
        } else if(energy > ParticleInfo(PID::pionp()).Mass()) {
            if(IsSinglePion(mode, kinetic, rng.uniform(0.0, 1.0))) {
                outgoing = SinglePion(mode, rng, part1CM, part2CM);
            } else  if (energy > 2*ParticleInfo(PID::pionp()).Mass()) {
                outgoing = DoublePion(mode, rng, part1CM, part2CM);
            } else {
                outgoing = Elastic(mode, rng, kinetic, part1CM, part2CM);
            }
        } else 
            outgoing = Elastic(mode, rng, kinetic, part1CM, part2CM);
    } else 
        outgoing = Elastic(mode, rng, kinetic, part1CM, part2CM);

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
    // Rambo: 4*n = 16
    // Channel select: 1
    // Position generation: 2
    // Total Random numbers = 19
    constexpr size_t ndims = 19;
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
    return GetProbe(mode.first) -> IsCEX(mode.second, kinetic, ran);
}

bool nuchic::Metropolis::IsSinglePion(InteractionType mode, const double &kinetic, const double &ran) const {
    return GetProbe(mode.first) -> IsSinglePion(mode.second, kinetic, ran);
}
