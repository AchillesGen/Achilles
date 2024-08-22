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
        algorithm = [&](Cascade *cascade, size_t idx, Particles &particles) -> size_t {
            return cascade->BaseAlgorithm(idx, particles);
        };
        break;
    case Algorithm::MFP:
        algorithm = [&](Cascade *cascade, size_t idx, Particles &particles) -> size_t {
            return cascade->MFPAlgorithm(idx, particles);
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

    spdlog::debug("Cascade initialized with distance: {}, medium: {}, potential_prop: {}, "
                  "probability: {}, algorithm: {}",
                  distance, ToString(m_medium), m_potential_prop, m_probability_name,
                  ToString(alg));

    kickedIdxs.clear();
}

void Cascade::Kick(Event &event, const FourVector &energyTransfer,
                   const std::array<double, 2> &sigma) {
    std::vector<std::size_t> indices;

    // Interact with protons or neutrons according to their total cross section
    auto ddSigma = {sigma[0], sigma[1]};
    auto index = Random::Instance().SelectIndex(ddSigma);

    auto interactPID = index == 0 ? PID::proton() : PID::neutron();

    // Restrict to particles of chosen species
    for(std::size_t i = 0; i < event.Hadrons().size(); ++i) {
        if(event.Hadrons()[i].ID() == interactPID) indices.push_back(i);
    }

    // Kick a single particle from the list
    kickedIdxs.insert(Random::Instance().Pick(indices));
    auto kicked = &event.Hadrons()[*kickedIdxs.begin()];
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

    auto mom = m_nucleus->GenerateMomentum(position);
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
            fact = m_nucleus->GetPotential()->InMediumCorrectionNonRel(p1, p2, mass, position);

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
            fact = m_nucleus->GetPotential()->InMediumCorrectionNonRel(p1, p2, mass, position);

        xsecDiff = GetXSec(kickedPart, particles[idxDiff]) * fact;
    }

    double rhoSame = 0.0;
    double rhoDiff = 0.0;
    if(position < m_nucleus->Radius()) {
        // TODO: Adjust below to handle non-isosymmetric nuclei
        rhoSame = m_nucleus->Rho(position) * 2 * static_cast<double>(index_same.size()) /
                  static_cast<double>(particles.size());
        rhoDiff = m_nucleus->Rho(position) * 2 * static_cast<double>(index_diff.size()) /
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
    kickedIdxs.clear();
    integrators.clear();
}

std::set<size_t> Cascade::InitializeIntegrator(Event &event) {
    std::set<size_t> notCaptured{};
    for(auto idx : kickedIdxs) {
        if(m_potential_prop && m_nucleus->GetPotential()->Hamiltonian(
                                   event.Hadrons()[idx].Momentum().P(),
                                   event.Hadrons()[idx].Position().P()) < Constant::mN) {
            event.Hadrons()[idx].Status() = ParticleStatus::captured;
        } else {
            AddIntegrator(idx, event.Hadrons()[idx]);
            notCaptured.insert(idx);
        }
    }
    kickedIdxs = notCaptured;
    return kickedIdxs;
}

void Cascade::UpdateKicked(Particles &particles, std::set<size_t> &newKicked) {
    for(size_t idx = 0; idx < particles.size(); ++idx) {
        if(particles[idx].Status() == ParticleStatus::propagating) {
            if(m_potential_prop &&
               m_nucleus->GetPotential()->Hamiltonian(
                   particles[idx].Momentum().P(), particles[idx].Position().P()) < Constant::mN) {
                particles[idx].Status() = ParticleStatus::captured;
            } else {
                AddIntegrator(idx, particles[idx]);
                newKicked.insert(idx);
            }
        }
    }
}

void Cascade::Validate(const Particles &particles) {
    for(auto particle : particles) {
        if(particle.Status() == ParticleStatus::propagating) {
            for(auto p : particles) spdlog::error("{}", p);
            throw std::runtime_error("Cascade has failed. Insufficient max steps.");
        }
    }
}

void Cascade::Evolve(achilles::Event &event, Nucleus *nucleus, const std::size_t &maxSteps) {
    // Set all propagating particles as kicked for the cascade
    for(size_t idx = 0; idx < event.Hadrons().size(); ++idx) {
        if(event.Hadrons()[idx].Status() == ParticleStatus::propagating) SetKicked(idx);
        else if(event.Hadrons()[idx].Status() == ParticleStatus::background) {
            auto mom3 = ThreeVector(nucleus->GenerateMomentum(event.Hadrons()[idx].Position().Magnitude()));
            auto mass = event.Hadrons()[idx].Info().Mass();
            auto energy = sqrt(mom3*mom3 + mass*mass);
            event.Hadrons()[idx].Momentum() = {mom3,energy};
        }
    }

    // Run the cascade
    m_nucleus = nucleus;
    Particles &particles = event.Hadrons();
    kickedIdxs = InitializeIntegrator(event);

    for(std::size_t step = 0; step < maxSteps; ++step) {
        // Stop loop if no particles are propagating
        if(kickedIdxs.size() == 0) break;

        // Adapt time step
        AdaptiveStep(particles, distance);

        // Make local copy of kickedIdxs
        std::set<size_t> newKicked{};
        for(auto idx : kickedIdxs) {
            Particle *kickNuc = &particles[idx];
            spdlog::trace("Kicked ID: {}, Particle: {}", idx, *kickNuc);

            // Update formation zones
            if(kickNuc->InFormationZone()) {
                Propagate(idx, kickNuc, distance);
                kickNuc->UpdateFormationZone(timeStep);
                newKicked.insert(idx);
                continue;
            }

            auto hitIdx = algorithm(this, idx, particles);
            if(hitIdx == SIZE_MAX) {
                newKicked.insert(idx);
                continue;
            }
            FinalizeMomentum(event, particles, idx, hitIdx);
        }
        // Updated the kicked list
        // NOTE: Needs to be here in case there are two interactions in the same time step
        UpdateKicked(particles, newKicked);

        // Replace kicked indices with new list
        kickedIdxs = newKicked;

        // After step checks
        Escaped(particles);
    }

    Validate(particles);
    Reset();
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
        SymplecticIntegrator(part.Position(), part.Momentum().Vec3(), m_nucleus->GetPotential(),
                             dHamiltonian_dr, dHamiltonian_dp, omega);
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

// TODO: Switch to using reference_wrapper to reduce copying
void Cascade::MeanFreePath(Event &event, Nucleus *nucleus, const std::size_t &maxSteps) {
    m_nucleus = nucleus;
    Particles particles = event.Hadrons();

    if(kickedIdxs.size() != 1) {
        throw std::runtime_error("MeanFreePath: only one particle should be kicked.");
    }

    auto idx = *kickedIdxs.begin();
    Particle *kickNuc = &particles[idx];

    // Initialize symplectic integrator
    AddIntegrator(idx, particles[idx]);
    if(m_potential_prop && m_nucleus->GetPotential()->Hamiltonian(
                               kickNuc->Momentum().P(), kickNuc->Position().P()) < Constant::mN) {
        kickNuc->Status() = ParticleStatus::captured;
        event.Hadrons() = particles;
        Reset();
        return;
    }

    if((kickNuc->Status() != ParticleStatus::internal_test) &&
       (kickNuc->Status() != ParticleStatus::external_test)) {
        throw std::runtime_error("MeanFreePath: kickNuc must have status -3 "
                                 "in order to accumulate DistanceTraveled.");
    }
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
        if(kickNuc->Position().Magnitude() >= nucleus->Radius() &&
           kickNuc->Status() == ParticleStatus::external_test) {
            if(kickNuc->Position().Pz() > nucleus->Radius()) break;

        } else if(kickNuc->Position().Magnitude() >= nucleus->Radius()) {
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

        // Did we *really* hit? Finalize momentum, check for Pauli blocking.
        FinalizeMomentum(event, particles, idx, hitIdx);
        // Stop as soon as we hit anything
        if(particles[idx].Status() == ParticleStatus::interacted) break;
    }
    event.Hadrons() = particles;
    Reset();
}

// TODO: Switch to using reference_wrapper to reduce copying
void Cascade::MeanFreePath_NuWro(Event &event, Nucleus *nucleus, const std::size_t &maxSteps) {
    m_nucleus = nucleus;
    Particles particles = event.Hadrons();

    if(kickedIdxs.size() != 1) {
        std::runtime_error("MeanFreePath: only one particle should be kicked.");
    }

    auto idx = *kickedIdxs.begin();
    Particle *kickNuc = &particles[idx];

    // Initialize symplectic integrator
    AddIntegrator(idx, particles[idx]);
    if(m_potential_prop && m_nucleus->GetPotential()->Hamiltonian(
                               kickNuc->Momentum().P(), kickNuc->Position().P()) < Constant::mN) {
        kickNuc->Status() = ParticleStatus::captured;
        event.Hadrons() = particles;
        Reset();
        return;
    }

    if(kickNuc->Status() != ParticleStatus::internal_test) {
        std::runtime_error("MeanFreePath: kickNuc must have status -3 "
                           "in order to accumulate DistanceTraveled.");
    }

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
        // Did we *really* hit? Finalize momentum, check for Pauli blocking.
        FinalizeMomentum(event, particles, idx, hitIdx);
        // Stop as soon as we hit anything
        if(particles[idx].Status() == ParticleStatus::interacted) break;
    }

    event.Hadrons() = particles;
    Reset();
}

// TODO: Rewrite to have the logic built into the Nucleus class
void Cascade::Escaped(Particles &particles) {
    const auto radius = m_nucleus->Radius();
    for(auto it = kickedIdxs.begin(); it != kickedIdxs.end();) {
        // Nucleon outside nucleus (will check if captured or escaped after
        // cascade)
        auto particle = &particles[*it];
        spdlog::debug("Kicked particle: {} at idx = {}", *particle, *it);
        if(particle->Status() != ParticleStatus::propagating &&
           particle->Status() != ParticleStatus::external_test) {
            spdlog::info("Particle: {}", *particle);
            throw std::domain_error("Invalid Particle in kicked list");
        }
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
            // TODO: Figure out how to appropriately handle escaping vs. capturing
            // It should not be returned to the background, since this can lead to
            // interactions with other particles again
            if(energy > 0 || !particle->Info().IsNucleon())
                particle->Status() = ParticleStatus::final_state;
            else
                particle->Status() = ParticleStatus::captured;
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

/// Determine whether "test_point" is between two parallel planes defined by
/// the normal vector "plane_normal". With the point "on_plane" on the closer plane
/// and dist the distance from the start plane to the end plane.
bool Cascade::BetweenPlanes(const ThreeVector &test_point, const ThreeVector &on_plane,
                            const ThreeVector &plane_normal, double dist) const noexcept {
    // Get signed distance to plane
    auto signed_dist = plane_normal.Dot(test_point - on_plane);

    // Determine if point is between planes
    return signed_dist > 0 && signed_dist < dist;
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
    auto distance2 = (point2 - point1).Dot(normedMomentum);

    // Build results vector
    for(std::size_t i = 0; i < particles.size(); ++i) {
        // TODO: Should particles propagating be able to interact with
        //       other propagating particles?
        if(particles[i].Status() != ParticleStatus::background) continue;
        // if(i == idx) continue;
        //  if(particles[i].InFormationZone()) continue;
        if(!BetweenPlanes(particles[i].Position(), point1, normedMomentum, distance2)) continue;
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
    // TODO: Clean this up so we don't have to calculate the cross section multiple times
    auto fact = InMediumCorrection(particle1, particle2);
    return m_interactions.TotalCrossSection(particle1, particle2) * fact;
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

void Cascade::FinalizeMomentum(Event &event, Particles &particles, size_t idx1,
                               size_t idx2) noexcept {
    Particle &particle1 = particles[idx1];
    Particle &particle2 = particles[idx2];

    if(Absorption(event, particle1, particle2)) return;

    // NOTE: No need to deal with in-medium effects since they are just an overall scaling
    auto modes = m_interactions.CrossSection(particle1, particle2);
    auto mode = m_interactions.SelectChannel(modes, Random::Instance().Uniform(0.0, 1.0));

    auto particles_out =
        m_interactions.GenerateMomentum(particle1, particle2, mode.particles, Random::Instance());

    spdlog::trace("Outgoing particles:");
    for(const auto &part : particles_out) spdlog::trace("- {}", part);

    // Check for Pauli Blocking
    bool hit = true;
    for(const auto &part : particles_out) hit &= !PauliBlocking(part);

    if(hit) {

        //spdlog::info("fermigas weight = {}", m_nucleus->FermiGasWeight(particle2));
        event.Weight() *= m_nucleus->FermiGasWeight(particle2);

        Particles initial_part, final_part;
        particle1.Status() = ParticleStatus::interacted;
        particle2.Status() = ParticleStatus::interacted;
        initial_part.push_back(particle1);
        initial_part.push_back(particle2);

        // Ensure outgoing particles are propagating and add to list of particles in event
        // and assign formation zone
        for(auto &part : particles_out) {
            part.Status() = ParticleStatus::propagating;
            part.SetFormationZone(particle1.Momentum(), part.Momentum());
            particles.push_back(part);
            final_part.push_back(particles.back());

            // TODO: Work out the queue to minimize number of calls to AllowedInteractions
            // Update queue of closest approach times
            // for(size_t i = 0; i < particles.size(); ++i) {
            //     if(particles[i].Status() != ParticleStatus::background) continue;
            //     double closest = ClosestApproach(part, particles[i]);
            //     spdlog::trace("Closest approach time({}) = {}", i, closest+currentTime);
            //     if(closest > 0) { m_time_steps.push({closest+currentTime, {particles.size()-1,
            //     i}}); }
            // }
        }

        // Add interaction to the event history
        // What do we use for the position? (How about average positions?)
        auto average_position = (particle1.Position() + particle2.Position()) / 2.0;
        event.History().AddVertex(average_position, initial_part, final_part,
                                  EventHistory::StatusCode::cascade);
    }
}

// TODO: Rewrite to have most of the logic built into the Nucleus class?
bool Cascade::PauliBlocking(const Particle &particle) const noexcept {
    if(particle.ID() != PID::proton() && particle.ID() != PID::neutron()) return false;
    double position = particle.Position().Magnitude();
    return particle.Momentum().Vec3().Magnitude() < m_nucleus->FermiMomentum(position);
}

double Cascade::InMediumCorrection(const Particle &particle1, const Particle &particle2) const {
    if(m_medium != InMedium::NonRelativistic) return 1;

    auto p1 = particle1.Momentum();
    auto p2 = particle2.Momentum();
    double mass = particle1.Info().Mass();
    auto pos_p1 = particle1.Position();
    auto pos_p2 = particle2.Position();
    double position1 = pos_p1.Magnitude();
    double position2 = pos_p2.Magnitude();
    double position3 = (pos_p1 + pos_p2).Magnitude();
    return m_nucleus->GetPotential()->InMediumCorrectionNonRel(p1, p2, mass, position1, position2,
                                                               position3);
}

bool Cascade::Absorption(Event &event, Particle &particle1, Particle &particle2) noexcept {
    // Check if absorption is possible
    // NOTE: Particle1 is always a propagating particle
    // NOTE: Particle2 is always a background particle
    if(!particle1.Info().IsPion() || !particle2.Info().IsNucleon()) return false;

    Particles &particles = event.Hadrons();

    Particle initial_pion;
    Particle incoming_nucleon1;
    Particle final_nucleon1;
    Particle final_nucleon2;

    int initial_charge = 0;
    int final_charge = 0;

    // Look for a previous pi-N Elastic or CE scattering
    auto node = event.History().FindNodeOut(particle1);
    if(node == nullptr || !node->IsCascade() || node->ParticlesOut().size() != 2) return false;

    bool CE = false;
    bool elastic = false;

    bool init_pion = false;
    bool init_nuc = false;
    bool out_nuc = false;
    // These are in initial particles of the previous pi-N scattering
    for(auto &in_part : node->ParticlesIn()) {
        if(in_part.Info().IsNucleon()) {
            // Found initial nucleon1
            incoming_nucleon1 = in_part;
            init_nuc = true;
        } else if(in_part.Info().IsPion()) {
            // Found initial pion
            initial_pion = in_part;
            init_pion = true;
        }
    }

    // These are the final particles of the previous pi-N scattering
    // This pion is the intermediate
    // This nucleon is the final nucleon1
    for(auto &out_part : node->ParticlesOut()) {
        if(out_part.Info().IsNucleon()) {
            // Found final nucleon 1
            final_nucleon1 = out_part;
            out_nuc = true;
            if(final_nucleon1.ID() == incoming_nucleon1.ID())
                elastic = true;
            else
                CE = true;
        }
    }

    // Makes sure that the we had pi-N -> pi-N
    if(!(init_pion && init_nuc && out_nuc)) return false;

    auto abs_prob = 1.0;
    if(Random::Instance().Uniform(0.0, 1.0) > abs_prob) return false;

    if(CE) {
        // We charge exchanged the 1st interaction
        final_nucleon2 = particle2;
    }

    if(elastic) {
        // We need to charge exchange the 2nd interaction
        if(particle2.ID() == PID::proton())
            final_nucleon2 = {PID::neutron(), particle2.Momentum(), particle2.Position()};
        if(particle2.ID() == PID::neutron())
            final_nucleon2 = {PID::proton(), particle2.Momentum(), particle2.Position()};
    }

    // Add up initial and final charges
    initial_charge = initial_pion.Info().IntCharge() + incoming_nucleon1.Info().IntCharge() +
                     particle2.Info().IntCharge();
    final_charge = final_nucleon1.Info().IntCharge() + final_nucleon2.Info().IntCharge();

    if(initial_charge != final_charge) return false;

    // Let's absorb the pion!
    // First let's get the momentum of the incoming pion for the current interaction
    auto combined_momentum = particle1.Momentum().Vec3() + particle2.Momentum().Vec3();

    // Now let's reset the past node and fill it with our initial system positions
    auto new_average_position =
        (initial_pion.Position() + incoming_nucleon1.Position() + final_nucleon2.Position()) / 3.0;

    final_nucleon2.Momentum() = {
        combined_momentum, sqrt(pow(final_nucleon2.Mass(), 2) + combined_momentum.Magnitude2())};

    final_nucleon2.Status() = ParticleStatus::propagating;

    // TODO : Do we need pauli blocking for absorption?
    // Check for Pauli Blocking
    // if(!PauliBlocking(final_nucleon2) || !PauliBlocking(final_nucleon1)) return false;

    particles.push_back(final_nucleon2);
    node->ResetParticles();
    node->SetPosition(new_average_position);
    node->AddIncoming(initial_pion);
    node->AddIncoming(incoming_nucleon1);
    node->AddIncoming(particle2);
    node->AddOutgoing(final_nucleon1);
    node->AddOutgoing(particles.back());

    particle1.Status() = ParticleStatus::interacted;
    particle2.Status() = ParticleStatus::interacted;

    return true;
}
