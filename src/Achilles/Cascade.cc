#include <map>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"

#include "Achilles/Cascade.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/Exception.hh"
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
        algorithm = [&](Cascade *cascade, size_t idx, Event &event) -> size_t {
            return cascade->BaseAlgorithm(idx, event);
        };
        break;
    case Algorithm::MFP:
        algorithm = [&](Cascade *cascade, size_t idx, Event &event) -> size_t {
            return cascade->MFPAlgorithm(idx, event);
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

std::size_t Cascade::GetInter(Particles &, const Particle &, double &) {
    throw std::logic_error("Cascade: MFPAlgorithm is currently not working with pions");
    /*
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
    */
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
        if(particles[idx].IsPropagating()) {
            if(m_potential_prop &&
               m_nucleus->GetPotential()->Hamiltonian(
                   particles[idx].Momentum().P(), particles[idx].Position().P()) < Constant::mN) {
                particles[idx].Status() = ParticleStatus::captured;
            } else {
                spdlog::trace("Adding particle: {}", particles[idx]);
                AddIntegrator(idx, particles[idx]);
                newKicked.insert(idx);
            }
        }
    }
}

void Cascade::Validate(Event &event) {
    spdlog::trace("End of Event Status:");
    for(size_t idx = 0; idx < event.Hadrons().size(); ++idx) {
        auto &particle = event.Hadrons()[idx];
        if(particle.Status() == ParticleStatus::propagating) {
            for(auto p : event.Hadrons()) spdlog::error("{}", p);
            achilles::PrintVisitor visitor;
            event.History().WalkHistory(visitor);
            spdlog::error("Event history: {}", visitor.data);
            throw AchillesCascadeError("Cascade has failed. Insufficient max steps.");
        }
        if(particle.Info().IsResonance() && particle.Status() == ParticleStatus::final_state) {
            static constexpr size_t nwarns = 10;
            static size_t iwarns = 0;
            if(iwarns < nwarns) {
                spdlog::warn("Cascade: Resonance did not decay before escaping, decaying now");
            } else if(iwarns == nwarns) {
                spdlog::warn(
                    "Cascade: Reached maximum escaped resonance warnings, suppressing the rest");
            }
            iwarns++;

            // Force decay
            auto particles_out = m_decays.Decay(particle);
            event.Hadrons()[idx].Status() = ParticleStatus::decayed;

            // Ensure outgoing particles are propagating and add to list of particles in event
            // and assign formation zone
            std::vector<Particle> final;
            for(auto &out : particles_out) {
                out.Status() = ParticleStatus::final_state;
                out.Position() = event.Hadrons()[idx].Position();
                event.Hadrons().push_back(out);
                final.push_back(event.Hadrons().back());
            }

            // Add decay to the event history
            event.History().AddVertex(event.Hadrons()[idx].Position(), {event.Hadrons()[idx]},
                                      final, EventHistory::StatusCode::decay);
        }
    }
}

void Cascade::PropagateAll(Particles &particles, double step) const {
    for(auto &particle : particles) {
        if(particle.IsPropagating()) { particle.Propagate(step); }
    }
}

void Cascade::Evolve(achilles::Event &event, Nucleus *nucleus,
                     [[maybe_unused]] const std::size_t &maxSteps) {
    // Set all propagating particles as kicked for the cascade
    for(size_t idx = 0; idx < event.Hadrons().size(); ++idx) {
        if(event.Hadrons()[idx].Status() == ParticleStatus::propagating) SetKicked(idx);
    }

    // Run the cascade
    currentTime = 0;
    m_nucleus = nucleus;
    Particles &particles = event.Hadrons();
    kickedIdxs = InitializeIntegrator(event);
    for(const auto &kicked : kickedIdxs) {
        for(size_t i = 0; i < particles.size(); ++i) {
            if(particles[i].Status() != ParticleStatus::background) continue;
            double closest = ClosestApproach(particles[kicked], particles[i]);
            spdlog::debug("Closest approach time({}, {}) = {}", kicked, i, closest);
            if(closest > 0) { m_time_steps.push({closest, {kicked, i}}); }
        }
    }

    // while(!m_time_steps.empty()) {
    //     timeStep = m_time_steps.top().time - currentTime;
    //     auto [kicked, hit] = m_time_steps.top().idxs;
    //     m_time_steps.pop();
    //     PropagateAll(particles, timeStep);
    //     currentTime += timeStep;
    //     if(HasInteraction(particles, kicked, hit)) {
    //         spdlog::trace("Kicked = {}", kicked);
    //         spdlog::trace("Hit = {}", hit);
    //         FinalizeMomentum(event, particles, kicked, hit);
    //         // tmp = m_time_steps;
    //         // pop_println("After Interaction", tmp);
    //         // if(m_time_steps.empty()) throw;
    //     }
    //     Escaped(particles);
    // }

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

            // Attempt decay of unstable particles
            if(!kickNuc->Info().IsStable()) {
                if(Decay(event, idx)) continue;
            }

            // Update formation zones
            if(kickNuc->InFormationZone() && !kickNuc->Info().IsPion()) {
                Propagate(idx, kickNuc);
                continue;
            }

            auto hitIdx = algorithm(this, idx, event);
            if(hitIdx == SIZE_MAX) { continue; }
            spdlog::trace("Kicked = {}", idx);
            spdlog::trace("Hit = {}", hitIdx);
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

    Validate(event);
    Reset();
}

size_t Cascade::BaseAlgorithm(size_t idx, Event &event) {
    // Get allowed interactions
    auto dist2 = AllowedInteractions(event.Hadrons(), idx);
    if(dist2.size() == 0) { return SIZE_MAX; }

    // Get interaction
    return Interacted(event, idx, dist2);
}

size_t Cascade::MFPAlgorithm(size_t idx, Event &event) {
    double step_prop = distance;
    Particle *kickNuc = &event.Hadrons()[idx];
    auto hitIdx = GetInter(event.Hadrons(), *kickNuc, step_prop);
    PropagateSpace(idx, kickNuc, step_prop);
    return hitIdx;
}

void Cascade::AddIntegrator(size_t idx, const Particle &part) {
    // Only add if not already there
    if(integrators.find(idx) != integrators.end()) return;

    // Otherwise setup integrator
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

void Cascade::Propagate(size_t idx, Particle *kickNuc) {
    if(m_potential_prop) {
        integrators[idx].Step<2>(timeStep);
        double energy = sqrt(pow(kickNuc->Info().Mass(), 2) + integrators[idx].P().P2());
        FourVector mom{integrators[idx].P(), energy};
        kickNuc->SetMomentum(mom);
        auto pos_old = kickNuc->Position();
        kickNuc->SetPosition(integrators[idx].Q());
        auto pos_new = kickNuc->Position();
        kickNuc->DistanceTraveled() += (pos_new - pos_old).Magnitude();
    } else {
        kickNuc->Propagate(timeStep);
    }
}

void Cascade::PropagateSpace(size_t idx, Particle *kickNuc, double step) {
    auto beta = kickNuc->Beta().Magnitude();
    timeStep = step / beta;
    Propagate(idx, kickNuc);
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
            Propagate(idx, kickNuc);
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
        auto hitIdx = Interacted(event, idx, nearby_particles);
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
        PropagateSpace(idx, kickNuc, step_prop);
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
    for(auto &particle : particles) {
        if(!particle.IsPropagating()) continue;
        if(particle.Status() == ParticleStatus::external_test) {
            if(particle.Position().Z() < radius) continue;
        }
        // TODO: Use the code from src/Achilles/Nucleus.cc:108 to properly
        // handle
        //       escape vs. capture and mometum changes
        constexpr double potential = 10.0;
        const double energy = particle.Momentum().E() - Constant::mN - potential;
        if(particle.Position().Magnitude2() > pow(radius, 2)) {
            spdlog::debug("Particle: {} escaping", particle);
            // TODO: Figure out how to appropriately handle escaping vs. capturing
            // It should not be returned to the background, since this can lead to
            // interactions with other particles again
            if(energy > 0 || !particle.Info().IsNucleon())
                particle.Status() = ParticleStatus::final_state;
            else
                particle.Status() = ParticleStatus::captured;
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

    timeStep = stepDistance / beta;
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
    Propagate(idx, &particles[idx]);
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

bool Cascade::HasInteraction(Event &event, size_t kicked, size_t hit) const {
    auto &particles = event.Hadrons();
    if(!particles[kicked].IsPropagating() ||
       particles[hit].Status() != ParticleStatus::background || particles[kicked].InFormationZone())
        return false;
    // Get distance between points
    const auto distance2 = (particles[kicked].Position() - particles[hit].Position()).Magnitude2();
    const auto xsec = GetXSec(event, kicked, hit);
    spdlog::trace("dist2[{}, {}] = {}", kicked, hit, distance2);
    const double prob = probability(distance2, xsec / 10);
    if(Random::Instance().Uniform(0.0, 1.0) < prob) {
        spdlog::debug("Accepted interaction with: {}, {} => distance2: {}, xsec: {}, prob: {}",
                      kicked + 1, hit + 1, distance2, xsec, prob);
        return true;
    }
    return false;
}

double Cascade::GetXSec(Event &event, size_t idx1, size_t idx2) const {
    // TODO: Clean this up so we don't have to calculate the cross section multiple times
    // TODO: Handle the in-medium effects for different processes
    auto fact = InMediumCorrection(event.Hadrons()[idx1], event.Hadrons()[idx2]);
    return m_interactions.TotalCrossSection(event, idx1, idx2) * fact;
}

/// Decide whether or not an interaction occured.
/// The total probability is normalized to the cross section "sigma" when
/// integrated over the plane.
std::size_t Cascade::Interacted(Event &event, size_t kicked,
                                const InteractionDistances &dists) noexcept {
    for(auto dist : dists) {
        // Cross section in mb
        const double xsec = GetXSec(event, kicked, dist.first);
        // 1 barn = 100 fm^2, so 1 mb = 0.1 fm^2.
        // Thus: (xsec [mb]) x (0.1 [fm^2]/ 1 [mb]) = 0.1 xsec [fm^2]
        // dist.second is fm^2; factor of 10 converts mb to fm^2
        const double prob = probability(dist.second, xsec / 10);
        spdlog::trace("dist2[{}, {}] = {}", kicked, dist.first, dist.second);
        if(Random::Instance().Uniform(0.0, 1.0) < prob) return dist.first;
    }

    return SIZE_MAX;
}

void Cascade::FinalizeMomentum(Event &event, Particles &particles, size_t idx1,
                               size_t idx2) noexcept {
    Particle &particle1 = particles[idx1];
    Particle &particle2 = particles[idx2];

    // NOTE: No need to deal with in-medium effects since they are just an overall scaling
    auto modes = m_interactions.CrossSection(event, idx1, idx2);
    auto mode = m_interactions.SelectChannel(modes, Random::Instance().Uniform(0.0, 1.0));

    spdlog::debug("Selected Mode: {}, {} -> {}", particle1.ID(), particle2.ID(),
                  fmt::join(mode.particles, ","));

    auto particles_out =
        m_interactions.GenerateMomentum(particle1, particle2, mode.particles, Random::Instance());

    spdlog::debug("Outgoing particles:");
    for(const auto &part : particles_out) spdlog::debug("- {}", part);

    // Check for Pauli Blocking
    bool hit = true;
    // bool pionIS = false;
    // bool baryonFS = true;

    // Did we start with pions
    // if(particle1.Info().IsPion() || particle2.Info().IsPion()) pionIS = true;

    for(const auto &part : particles_out) {
        hit &= !PauliBlocking(part);
        // Are there any mesons in the final state?
        // if(!part.Info().IsBaryon()) baryonFS = false;
    }

    // If there are pions in in the initial state
    // and no mesons in the final state
    /*if(pionIS && baryonFS) {
        hit = true;
        spdlog::debug("Checking pion abs mom");
        spdlog::debug("Particle1 ({}): {}, Particle2 ({}): {}", particle1.ID(),
                      particle1.Momentum().Vec3().Magnitude(), particle2.ID(),
                      particle2.Momentum().Vec3().Magnitude());
        spdlog::debug("Pion abs FS mom: {}, {}", particles_out[0].Momentum().Vec3().Magnitude(),
                      particles_out[1].Momentum().Vec3().Magnitude());
        if(PauliBlocking(particles_out[0]) || PauliBlocking(particles_out[1])) {
            spdlog::info("Pauliblocked in absorption");
            spdlog::info("outgoing part 1 ({}): {}",
    particles_out[0].ID(),particles_out[0].Momentum().Vec3().Magnitude()); spdlog::info("part 1
    Fermi momentum = {}",m_nucleus->FermiMomentum(particles_out[0].Position().Magnitude()));
            spdlog::info("outgoing part 2 ({}): {}",
    particles_out[1].ID(),particles_out[1].Momentum().Vec3().Magnitude()); spdlog::info("part 2
    Fermi momentum = {}",m_nucleus->FermiMomentum(particles_out[1].Position().Magnitude()));
        }
    }*/

    // for(auto &part : event.Hadrons()) {
    //     if(part.Status() == ParticleStatus::absorption_partner)
    //         part.Status() = hit ? ParticleStatus::interacted : ParticleStatus::background;
    // }

    if(hit) {
        Particles initial_part, final_part;
        particle1.Status() = ParticleStatus::interacted;
        particle2.Status() = ParticleStatus::interacted;
        initial_part.push_back(particle1);
        initial_part.push_back(particle2);

        // Ensure outgoing particles are propagating and add to list of particles in event
        // and assign formation zone
        for(auto &part : particles_out) {
            part.Status() = ParticleStatus::propagating;
            if(part.Info().IsNucleon())
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
        // TODO: What do we use for the position? (How about average positions?)
        // TODO: How to best include the absorp_partner
        auto average_position = (particle1.Position() + particle2.Position()) / 2.0;
        event.History().AddVertex(average_position, initial_part, final_part,
                                  EventHistory::StatusCode::cascade);
    }
}

// TODO: Rewrite to have most of the logic built into the Nucleus class?
bool Cascade::PauliBlocking(const Particle &particle) const noexcept {
    if(!particle.Info().IsNucleon()) return false;
    double position = particle.Position().Magnitude();
    return particle.Momentum().Vec3().Magnitude() < m_nucleus->FermiMomentum(position);
}

double Cascade::InMediumCorrection(const Particle &particle1, const Particle &particle2) const {
    if(m_medium != InMedium::NonRelativistic) return 1;

    if(!particle1.Info().IsNucleon() || !particle2.Info().IsNucleon()) return 1;

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

bool Cascade::Decay(Event &event, size_t idx) const {
    auto part = event.Hadrons()[idx];
    double lifetime = Constant::HBARC / part.Info().Width();
    double decay_prob = exp(-timeStep / lifetime);
    // Should we attempt a decay in this time step
    spdlog::trace("decay prob = {}, {}, {}", decay_prob, timeStep, lifetime);
    if(Random::Instance().Uniform(0.0, 1.0) < decay_prob) return false;

    // Look up in decay handler
    auto particles_out = m_decays.Decay(part);

    for(const auto &out : particles_out) {
        if(PauliBlocking(out)) return false;
    }

    event.Hadrons()[idx].Status() = ParticleStatus::decayed;

    // Ensure outgoing particles are propagating and add to list of particles in event
    // and assign formation zone
    std::vector<Particle> final;
    for(auto &out : particles_out) {
        if(std::isnan(out.Momentum()[0])) {
            spdlog::error("Nan momentum in decay");
            spdlog::error("Pin = {}, Pout = [{}, {}]", part, particles_out[0], particles_out[1]);
            throw AchillesCascadeError("Nan momentum in decay");
        }
        out.Status() = ParticleStatus::propagating;
        out.SetFormationZone(out.Momentum(), part.Momentum());
        event.Hadrons().push_back(out);
        final.push_back(event.Hadrons().back());
    }

    // Add decay to the event history
    event.History().AddVertex(part.Position(), {event.Hadrons()[idx]}, final,
                              EventHistory::StatusCode::decay);

    return true;
}
