#include <map>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"

#include "Achilles/Cascade.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Potential.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"

using namespace achilles;

Cascade::Cascade(InteractionHandler interactions, const ProbabilityType &prob, Algorithm alg,
                 const InMedium &medium, bool potential_prop, double dist)
    : distance(dist), m_interactions(std::move(interactions)), m_medium(medium),
      m_potential_prop(potential_prop) {
    switch(alg) {
    case Algorithm::Base:
        algorithm = [&](size_t idx, Particles &particles) -> size_t {
            return BaseAlgorithm(idx, particles);
        };
        break;
    case Algorithm::MFP:
        algorithm = [&](size_t idx, Particles &particles) -> size_t {
            return MFPAlgorithm(idx, particles);
        };
        break;
    }

    switch(prob) {
    case ProbabilityType::Gaussian:
        m_probability_name = "Gaussian";
        probability = [](const double &b2, const double &sigma) -> double {
            return exp(-M_PI * b2 / sigma);
        };
        break;
    case ProbabilityType::Pion:
        m_probability_name = "Pion";
        probability = [](const double &b2, const double &sigma) -> double {
            double b = sqrt(b2);
            // TODO: This does not work, we need to rethink this
            // return
            // (135_MeV*sigma)/Constant::HBARC/(2*M_PI*b)*exp(-135_MeV*b/Constant::HBARC);
            return exp(-sqrt(2 * M_PI / sigma) * b);
        };
        break;
    case ProbabilityType::Cylinder:
        m_probability_name = "Cylinder";
        probability = [](const double &b2, const double &sigma) -> double {
            return b2 < sigma / M_PI ? 1 : 0;
        };
    }

    kickedIdxs.resize(0);
}

void Cascade::Kick(std::shared_ptr<Nucleus> nucleus, const FourVector &energyTransfer,
                   const std::array<double, 2> &sigma) {
    std::vector<std::size_t> indices;

    // Interact with protons or neutrons according to their total cross section
    auto ddSigma = {sigma[0], sigma[1]};
    auto index = Random::Instance().SelectIndex(ddSigma);

    auto interactPID = index == 0 ? PID::proton() : PID::neutron();

    // Restrict to particles of chosen species
    for(std::size_t i = 0; i < nucleus->Nucleons().size(); ++i) {
        if(nucleus->Nucleons()[i].ID() == interactPID) indices.push_back(i);
    }

    // Kick a single particle from the list
    kickedIdxs.push_back(Random::Instance().Pick(indices));
    auto kicked = &nucleus->Nucleons()[kickedIdxs.back()];
    kicked->Status() = ParticleStatus::propagating;
    kicked->SetMomentum(kicked->Momentum() + energyTransfer);
}

std::size_t Cascade::GetInter(Particles &particles, const Particle &kickedPart,
                              double &stepDistance) {
    std::vector<std::size_t> index_same, index_diff;

    for(std::size_t i = 0; i < particles.size(); ++i) {
        if(particles[i].Status() != ParticleStatus::background) continue;
        if(particles[i].ID() == kickedPart.ID())
            index_same.push_back(i);
        else
            index_diff.push_back(i);
    }

    if(index_diff.size() == 0 && index_same.size() == 0) return SIZE_MAX;

    double position = kickedPart.Position().Magnitude();
    auto p1 = kickedPart.Momentum();
    double mass = kickedPart.Info().Mass();

    auto mom = localNucleus->GenerateMomentum(position);
    double energy = pow(kickedPart.Info().Mass(), 2);
    for(auto p : mom) energy += p * p;
    std::size_t idxSame = SIZE_MAX;
    double xsecSame = 0;
    if(index_same.size() != 0) {
        idxSame = Random::Instance().Pick(index_same);
        particles[idxSame].SetMomentum(FourVector(mom[0], mom[1], mom[2], sqrt(energy)));

        auto p2 = particles[idxSame].Momentum();
        double fact = 1.0;
        if(m_medium == InMedium::NonRelativistic)
            fact = localNucleus->GetPotential()->InMediumCorrectionNonRel(p1, p2, mass, position);

        xsecSame = GetXSec(kickedPart, particles[idxSame]) * fact;
    }

    auto otherMass = kickedPart.ID() == PID::proton() ? ParticleInfo(PID::neutron()).Mass()
                                                      : ParticleInfo(PID::proton()).Mass();
    energy = otherMass * otherMass;
    for(auto p : mom) energy += p * p;
    std::size_t idxDiff = SIZE_MAX;
    double xsecDiff = 0;
    if(index_diff.size() != 0) {
        idxDiff = Random::Instance().Pick(index_diff);
        particles[idxDiff].SetMomentum(FourVector(mom[0], mom[1], mom[2], sqrt(energy)));

        auto p2 = particles[idxSame].Momentum();
        double fact = 1.0;
        if(m_medium == InMedium::NonRelativistic)
            fact = localNucleus->GetPotential()->InMediumCorrectionNonRel(p1, p2, mass, position);

        xsecDiff = GetXSec(kickedPart, particles[idxDiff]) * fact;
    }

    double rhoSame = 0.0;
    double rhoDiff = 0.0;
    if(position < localNucleus->Radius()) {
        // TODO: Adjust below to handle non-isosymmetric nuclei
        rhoSame = localNucleus->Rho(position) * 2 * static_cast<double>(index_same.size()) /
                  static_cast<double>(particles.size());
        rhoDiff = localNucleus->Rho(position) * 2 * static_cast<double>(index_diff.size()) /
                  static_cast<double>(particles.size());
    }
    if(rhoSame <= 0.0 && rhoDiff <= 0.0) return SIZE_MAX;
    double lambda_tilde = 1.0 / (xsecSame / 10 * rhoSame + xsecDiff / 10 * rhoDiff);
    double lambda = -log(Random::Instance().Uniform(0.0, 1.0)) * lambda_tilde;

    if(lambda > stepDistance) return SIZE_MAX;

    stepDistance = lambda;
    double ichoice = Random::Instance().Uniform(0.0, 1.0);
    if(ichoice < xsecSame / (xsecSame + xsecDiff)) {
        particles[idxSame].SetPosition(kickedPart.Position());
        return idxSame;
    }

    particles[idxDiff].SetPosition(kickedPart.Position());
    return idxDiff;
}

void Cascade::Reset() {
    kickedIdxs.resize(0);
    integrators.clear();
}

std::vector<size_t> Cascade::InitializeIntegrator() {
    std::vector<size_t> notCaptured{};
    for(auto idx : kickedIdxs) {
        if(m_potential_prop && localNucleus->GetPotential()->Hamiltonian(
                                   localNucleus->Nucleons()[idx].Momentum().P(),
                                   localNucleus->Nucleons()[idx].Position().P()) < Constant::mN) {
            localNucleus->Nucleons()[idx].Status() = ParticleStatus::captured;
        } else {
            AddIntegrator(idx, localNucleus->Nucleons()[idx]);
            notCaptured.push_back(idx);
        }
    }
    kickedIdxs = notCaptured;
    return kickedIdxs;
}

void Cascade::UpdateKicked(bool hit, size_t idx, size_t hitIdx, Particle *kickNuc, Particle *hitNuc,
                           std::vector<size_t> &newKicked) {
    if(hit) {
        if(m_potential_prop &&
           localNucleus->GetPotential()->Hamiltonian(kickNuc->Momentum().P(),
                                                     kickNuc->Position().P()) < Constant::mN) {
            kickNuc->Status() = ParticleStatus::captured;
        } else {
            newKicked.push_back(idx);
        }
        if(m_potential_prop && localNucleus->GetPotential()->Hamiltonian(
                                   hitNuc->Momentum().P(), hitNuc->Position().P()) < Constant::mN) {
            hitNuc->Status() = ParticleStatus::captured;
        } else {
            newKicked.push_back(hitIdx);
            AddIntegrator(hitIdx, *hitNuc);
            hitNuc->Status() = ParticleStatus::propagating;
        }
    } else {
        newKicked.push_back(idx);
    }
}

void Cascade::Validate(const Particles &particles) {
    for(auto particle : particles) {
        if(particle.Status() == ParticleStatus::propagating) {
            for(auto p : particles) spdlog::error("{}", p);
            throw std::runtime_error("Cascade has failed. Insufficient max steps.");
        }
    }
    localNucleus->Nucleons() = particles;
}

void Cascade::Evolve(achilles::Event *event, const std::size_t &maxSteps) {
    // Set all propagating particles as kicked for the cascade
    for(size_t idx = 0; idx < event->Hadrons().size(); ++idx) {
        if(event->Hadrons()[idx].Status() == ParticleStatus::propagating) SetKicked(idx);
    }

    // Run the normal cascade
    Evolve(event->CurrentNucleus(), maxSteps);
}

size_t Cascade::BaseAlgorithm(size_t idx, Particles &particles) {
    // Get allowed interactions
    auto dist2 = AllowedInteractions(particles, idx);
    if(dist2.size() == 0) { return SIZE_MAX; }

    // Get interaction
    return Interacted(particles, particles[idx], dist2);
}

size_t Cascade::MFPAlgorithm(size_t idx, Particles &particles) {
    double step_prop = distance;
    Particle *kickNuc = &particles[idx];
    auto hitIdx = GetInter(particles, *kickNuc, step_prop);
    Propagate(idx, kickNuc, step_prop);
    return hitIdx;
}

void Cascade::Evolve(std::shared_ptr<Nucleus> nucleus, const std::size_t &maxSteps) {
    localNucleus = nucleus;
    Particles particles = nucleus->Nucleons();
    kickedIdxs = InitializeIntegrator();

    for(std::size_t step = 0; step < maxSteps; ++step) {
        // Stop loop if no particles are propagating
        if(kickedIdxs.size() == 0) break;

        // Adapt time step
        AdaptiveStep(particles, distance);

        // Make local copy of kickedIdxs
        std::vector<size_t> newKicked{};
        for(auto idx : kickedIdxs) {
            Particle *kickNuc = &particles[idx];
            spdlog::debug("Kicked ID: {}, Particle: {}", idx, *kickNuc);

            // Update formation zones
            if(kickNuc->InFormationZone()) {
                Propagate(idx, kickNuc, distance);
                kickNuc->UpdateFormationZone(timeStep);
                newKicked.push_back(idx);
                continue;
            }

            auto hitIdx = algorithm(idx, particles);
            if(hitIdx == SIZE_MAX) {
                newKicked.push_back(idx);
                continue;
            }
            Particle *hitNuc = &particles[hitIdx];
            bool hit = FinalizeMomentum(*kickNuc, *hitNuc);
            UpdateIntegrator(idx, kickNuc);

            UpdateKicked(hit, idx, hitIdx, kickNuc, hitNuc, newKicked);
            spdlog::debug("newKicked size = {}, {}", newKicked.size(), hit);
        }

        // Replace kicked indices with new list
        kickedIdxs = newKicked;

        // After step checks
        Escaped(particles);
    }

    Validate(particles);
    Reset();
}

void Cascade::AddIntegrator(size_t idx, const Particle &part) {
    static constexpr double omega = 20;
    auto dHamiltonian_dr = [&](const ThreeVector &q, const ThreeVector &p,
                               std::shared_ptr<Potential> potential) {
        auto vals = potential->operator()(p.P(), q.P());
        auto dpot_dp = potential->derivative_p(p.P(), q.P());

        auto mass_eff =
            achilles::Constant::mN + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
        double numerator = (vals.rscalar + achilles::Constant::mN) * dpot_dp.rscalar + p.P();
        double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
        return numerator / denominator * p / p.P() + dpot_dp.rvector * p / p.P();
    };
    auto dHamiltonian_dp = [&](const ThreeVector &q, const ThreeVector &p,
                               std::shared_ptr<Potential> potential) {
        auto vals = potential->operator()(p.P(), q.P());
        auto dpot_dp = potential->derivative_p(p.P(), q.P());

        auto mass_eff =
            achilles::Constant::mN + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
        double numerator = (vals.rscalar + achilles::Constant::mN) * dpot_dp.rscalar + p.P();
        double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
        return numerator / denominator * p / p.P() + dpot_dp.rvector * p / p.P();
    };
    integrators[idx] =
        SymplecticIntegrator(part.Position(), part.Momentum().Vec3(), localNucleus->GetPotential(),
                             dHamiltonian_dr, dHamiltonian_dp, omega);
}

void Cascade::UpdateIntegrator(size_t idx, Particle *kickNuc) {
    integrators[idx].State() = PSState(kickNuc->Position(), kickNuc->Momentum().Vec3());
}

void Cascade::Propagate(size_t idx, Particle *kickNuc, double step) {
    double localTimeStep = step / (kickNuc->Beta().Magnitude());
    if(m_potential_prop) {
        integrators[idx].Step<2>(localTimeStep);
        double energy = sqrt(pow(kickNuc->Info().Mass(), 2) + integrators[idx].P().P2());
        FourVector mom{integrators[idx].P(), energy};
        kickNuc->SetMomentum(mom);
        auto pos_old = kickNuc->Position();
        kickNuc->SetPosition(integrators[idx].Q());
        auto pos_new = kickNuc->Position();
        kickNuc->DistanceTraveled() += (pos_new - pos_old).Magnitude();
    } else {
        kickNuc->SpacePropagate(step);
    }
}

void Cascade::MeanFreePath(std::shared_ptr<Nucleus> nucleus, const std::size_t &maxSteps) {
    localNucleus = nucleus;
    Particles particles = nucleus->Nucleons();

    if(kickedIdxs.size() != 1) {
        throw std::runtime_error("MeanFreePath: only one particle should be kicked.");
    }

    auto idx = kickedIdxs[0];
    Particle *kickNuc = &particles[idx];

    // Initialize symplectic integrator
    AddIntegrator(idx, particles[idx]);
    if(m_potential_prop && localNucleus->GetPotential()->Hamiltonian(
                               kickNuc->Momentum().P(), kickNuc->Position().P()) < Constant::mN) {
        kickNuc->Status() = ParticleStatus::captured;
        nucleus->Nucleons() = particles;
        Reset();
        return;
    }

    if(kickNuc->Status() != ParticleStatus::internal_test) {
        throw std::runtime_error("MeanFreePath: kickNuc must have status -3 "
                                 "in order to accumulate DistanceTraveled.");
    }
    bool hit = false;
    for(std::size_t step = 0; step < maxSteps; ++step) {
        AdaptiveStep(particles, distance);

        // if(kickNuc->IsUnstable()) {
        // Check if decay
        // }

        if(kickNuc->InFormationZone()) {
            Propagate(idx, kickNuc, distance);
            kickNuc->UpdateFormationZone(timeStep);
            continue;
        }

        // Are we already outside nucleus?
        if(kickNuc->Position().Magnitude() >= nucleus->Radius()) {
            kickNuc->Status() = ParticleStatus::final_state;
            break;
        }
        AdaptiveStep(particles, distance);
        // Identify nearby particles which might interact
        auto nearby_particles = AllowedInteractions(particles, idx);
        if(nearby_particles.size() == 0) continue;
        // Did we hit?
        auto hitIdx = Interacted(particles, *kickNuc, nearby_particles);
        if(hitIdx == SIZE_MAX) continue;
        Particle *hitNuc = &particles[hitIdx];
        // Did we *really* hit? Finalize momentum, check for Pauli blocking.
        hit = FinalizeMomentum(*kickNuc, *hitNuc);
        // Stop as soon as we hit anything
        if(hit) break;
    }
    nucleus->Nucleons() = particles;
    Reset();
}

void Cascade::MeanFreePath_NuWro(std::shared_ptr<Nucleus> nucleus, const std::size_t &maxSteps) {
    localNucleus = nucleus;
    Particles particles = nucleus->Nucleons();

    if(kickedIdxs.size() != 1) {
        std::runtime_error("MeanFreePath: only one particle should be kicked.");
    }

    auto idx = kickedIdxs[0];
    Particle *kickNuc = &particles[idx];

    // Initialize symplectic integrator
    AddIntegrator(idx, particles[idx]);
    if(m_potential_prop && localNucleus->GetPotential()->Hamiltonian(
                               kickNuc->Momentum().P(), kickNuc->Position().P()) < Constant::mN) {
        kickNuc->Status() = ParticleStatus::captured;
        nucleus->Nucleons() = particles;
        Reset();
        return;
    }

    if(kickNuc->Status() != ParticleStatus::internal_test) {
        std::runtime_error("MeanFreePath: kickNuc must have status -3 "
                           "in order to accumulate DistanceTraveled.");
    }

    bool hit = false;
    for(std::size_t step = 0; step < maxSteps; ++step) {
        // Are we already outside nucleus?
        if(kickNuc->Position().Magnitude() >= nucleus->Radius()) {
            kickNuc->Status() = ParticleStatus::final_state;
            break;
        }
        // AdaptiveStep(particles, distance);
        double step_prop = distance;
        // Did we hit?
        auto hitIdx = GetInter(particles, *kickNuc, step_prop);
        Propagate(idx, kickNuc, step_prop);
        if(hitIdx == SIZE_MAX) continue;
        Particle *hitNuc = &particles[hitIdx];
        // Did we *really* hit? Finalize momentum, check for Pauli blocking.
        hit = FinalizeMomentum(*kickNuc, *hitNuc);
        // Stop as soon as we hit anything
        if(hit) break;
    }

    nucleus->Nucleons() = particles;
    Reset();
}

// TODO: Rewrite to have the logic built into the Nucleus class
void Cascade::Escaped(Particles &particles) {
    const auto radius = localNucleus->Radius();
    for(auto it = kickedIdxs.begin(); it != kickedIdxs.end();) {
        // Nucleon outside nucleus (will check if captured or escaped after
        // cascade)
        auto particle = &particles[*it];
        if(particle->Status() == ParticleStatus::background)
            throw std::domain_error("Invalid Particle in kicked list");
        // if(particle -> Status() == -2) {
        //     std::cout << particle -> Position().Pz() << " " << sqrt(radius2)
        //     << std::endl; std::cout << *particle << std::endl;
        // }
        // TODO: Use the code from src/Achilles/Nucleus.cc:108 to properly
        // handle
        //       escape vs. capture and mometum changes
        constexpr double potential = 10.0;
        const double energy = particle->Momentum().E() - Constant::mN - potential;
        if(particle->Position().Magnitude2() > pow(radius, 2) &&
           particle->Status() != ParticleStatus::external_test) {
            if(energy > 0)
                particle->Status() = ParticleStatus::final_state;
            else
                particle->Status() = ParticleStatus::background;
            it = kickedIdxs.erase(it);
        } else if(particle->Status() == ParticleStatus::external_test &&
                  particle->Position().Pz() > radius) {
            it = kickedIdxs.erase(it);
        } else {
            ++it;
        }
    }
}

/// Convert a time step in [fm] to [1/MeV].
/// timeStep = distance / max("betas of all kicked particles") / hbarc
void Cascade::AdaptiveStep(const Particles &particles, const double &stepDistance) noexcept {
    double beta = 0;
    for(auto idx : kickedIdxs) {
        if(particles[idx].Beta().Magnitude() > beta) beta = particles[idx].Beta().Magnitude();
    }

    timeStep = stepDistance / (beta * Constant::HBARC);
}

/// Determine whether "position" is between the two planes orthogonal to the
/// displacement "point2 - point1" and containing the points "point1" and
/// "point2,"" respectively.
bool Cascade::BetweenPlanes(const ThreeVector &position, const ThreeVector &point1,
                            const ThreeVector &point2) const noexcept {
    // Get distance between planes
    const ThreeVector dist = point2 - point1;

    // Determine if point is between planes
    return ((position - point1).Dot(dist) >= 0 && (position - point2).Dot(dist) <= 0);
}

/// Project the vector "position" onto the plane containing the point "planePt"
/// and orthogonal to the vector "planeVec" by subtracting off the projection
/// of "position - planePt" onto the normal vector "planeVec".
const ThreeVector Cascade::Project(const ThreeVector &position, const ThreeVector &planePt,
                                   const ThreeVector &planeVec) const noexcept {
    // Project point onto a plane
    const ThreeVector projection = (position - planePt).Dot(planeVec) * planeVec;
    return position - projection;
}

/// Get a sorted list of allowed InteractionDistances, i.e., of pairs
/// (index, distance). Allowed interactions are those between the planes
/// orthogonal to the momentum of the specified particle and located at
/// (i) the particle's current position and
/// (ii) the particle's position after propagating for "timeStep".
/// For example, in the diagram below, interaction with A is allowed, while
/// interactions with B and C are not allowed.
///    plane1  plane2
///      |   A   |
///      |       |
///      x---->  |
///      |       |   B
///  C   |       |
const InteractionDistances Cascade::AllowedInteractions(Particles &particles,
                                                        const std::size_t &idx) noexcept {
    InteractionDistances results;

    // Build planes
    const ThreeVector point1 = particles[idx].Position();
    Propagate(idx, &particles[idx], distance);
    const ThreeVector point2 = particles[idx].Position();
    auto normedMomentum = particles[idx].Momentum().Vec3().Unit();

    // Build results vector
    for(std::size_t i = 0; i < particles.size(); ++i) {
        // TODO: Should particles propagating be able to interact with
        //       other propagating particles?
        if(particles[i].Status() != ParticleStatus::background) continue;
        // if(i == idx) continue;
        //  if(particles[i].InFormationZone()) continue;
        if(!BetweenPlanes(particles[i].Position(), point1, point2)) continue;
        auto projectedPosition = Project(particles[i].Position(), point1, normedMomentum);
        // (Squared) distance in the direction orthogonal to the momentum
        double dist2 = (projectedPosition - point1).Magnitude2();

        results.push_back(std::make_pair(i, dist2));
    }

    // Sort array by distances
    std::sort(results.begin(), results.end(), sortPairSecond);

    return results;
}

double Cascade::GetXSec(const Particle &particle1, const Particle &particle2) const {
    auto p1 = particle1.Momentum();
    auto p2 = particle2.Momentum();
    double mass = particle1.Info().Mass();
    auto pos_p1 = particle1.Position();
    auto pos_p2 = particle2.Position();
    double position1 = pos_p1.Magnitude();
    double position2 = pos_p2.Magnitude();
    double position3 = (pos_p1 + pos_p2).Magnitude();
    double fact = 1.0;
    if(m_medium == InMedium::NonRelativistic)
        fact = localNucleus->GetPotential()->InMediumCorrectionNonRel(p1, p2, mass, position1,
                                                                      position2, position3);

    // TODO: Clean this up so we don't have to calculate the cross section multiple times
    auto modes = m_interactions.CrossSection(particle1, particle2);
    double xsec = 0;
    for(auto mode : modes) { xsec += mode.cross_section * fact; }
    return xsec;
}

/// Decide whether or not an interaction occured.
/// The total probability is normalized to the cross section "sigma" when
/// integrated over the plane.
std::size_t Cascade::Interacted(const Particles &particles, const Particle &kickedParticle,
                                const InteractionDistances &dists) noexcept {
    for(auto dist : dists) {
        // Cross section in mb
        const double xsec = GetXSec(kickedParticle, particles[dist.first]);
        // 1 barn = 100 fm^2, so 1 mb = 0.1 fm^2.
        // Thus: (xsec [mb]) x (0.1 [fm^2]/ 1 [mb]) = 0.1 xsec [fm^2]
        // dist.second is fm^2; factor of 10 converts mb to fm^2
        const double prob = probability(dist.second, xsec / 10);
        if(Random::Instance().Uniform(0.0, 1.0) < prob) return dist.first;
    }

    return SIZE_MAX;
}

bool Cascade::FinalizeMomentum(Particle &particle1, Particle &particle2) noexcept {
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    // TODO: Update this to the new interface
    auto pOut = m_interactions.GenerateMomentum(particle1, particle2,
                                                {particle1.ID(), particle2.ID()}, {0});

    // Assign momenta to particles
    particle1.SetMomentum(pOut[0].Momentum());
    particle2.SetMomentum(pOut[1].Momentum());

    // Check for Pauli Blocking
    bool hit = !(PauliBlocking(particle1) || PauliBlocking(particle2));

    if(hit) {
        // Assign formation zone
        particle1.SetFormationZone(p1Lab, pOut[0].Momentum());
        particle2.SetFormationZone(p2Lab, pOut[0].Momentum());

        // Hit nucleon is now propagating
        // Users are responsibile for updating the status externally as desired
    } else {
        // Assign momenta to particles
        particle1.SetMomentum(p1Lab);
        particle2.SetMomentum(p2Lab);
    }

    return hit;
}

// TODO: Rewrite to have most of the logic built into the Nucleus class?
bool Cascade::PauliBlocking(const Particle &particle) const noexcept {
    if(particle.ID() != PID::proton() || particle.ID() != PID::neutron()) return false;
    double position = particle.Position().Magnitude();
    return particle.Momentum().Vec3().Magnitude() < localNucleus->FermiMomentum(position);
}
