#include <random>
#include <iostream>

#include "spdlog/spdlog.h"

#include "nuchic/Cascade.hh"
#include "nuchic/Constants.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Utilities.hh"
#include "nuchic/Interactions.hh"

using namespace nuchic;

Particles Cascade::Kick(const Particles& particles, const FourVector& energyTransfer,
        const std::array<double, 2>& sigma) {
    std::vector<std::size_t> indices;

    auto ddSigma = {sigma[0], sigma[1]};
    std::size_t index = rng.variate<int, std::discrete_distribution>(ddSigma);

    int interactPID = index == 0 ? 2212 : 2112;

    for(std::size_t i = 0; i < particles.size(); ++i) {
        if(particles[i].PID() == interactPID) indices.push_back(i);
    }

    kickedIdxs.push_back(rng.pick(indices));
    Particles result = particles;
    result[kickedIdxs.back()].SetStatus(-1);
    result[kickedIdxs.back()].SetMomentum(result[kickedIdxs.back()].Momentum() + energyTransfer);

    return result;
}

void Cascade::Reset() {
    kickedIdxs.resize(0);
}

Particles Cascade::operator()(const Particles& _particles, const double& kf, const double& _radius2,
        const std::size_t& maxSteps) {

    Particles particles = _particles;
    fermiMomentum = kf;
    radius2= _radius2;
    for(std::size_t step = 0; step < maxSteps; ++step) {
        // Stop loop if no particles are propagating
        if(kickedIdxs.size() == 0) break;
        
        // Adapt time step
        AdaptiveStep(particles, distance);

        // Make local copy of kickedIdxs
        auto newKicked = kickedIdxs;
        for(auto idx : kickedIdxs) {
            Particle* kickNuc = &particles[idx];

            // Update formation zones
            if(kickNuc -> InFormationZone()) {
                kickNuc -> UpdateFormationZone(timeStep);
                kickNuc -> Propagate(timeStep);
                continue;
            }

            // Get allowed interactions
            auto dist2 = AllowedInteractions(particles, idx);
            if(dist2.size() == 0) continue;

            // Get interaction
            auto hitIdx = Interacted(particles, *kickNuc, dist2);
            if(hitIdx == -1) continue;
            Particle* hitNuc = &particles[hitIdx];

            // Finalize Momentum
            bool hit = FinalizeMomentum(*kickNuc, *hitNuc);
            if(hit) {
                newKicked.push_back(hitIdx);
                hitNuc->SetStatus(-1);
            }
        }

        // Replace kicked indices with new list
        kickedIdxs = newKicked;

        // After step checks
        //std::cout << "CURRENT STEP: " << step << std::endl;
        for(auto it = kickedIdxs.begin() ; it != kickedIdxs.end(); ) {
            // Nucleon outside nucleus (will check if captured or escaped after cascade)
            auto particle = &particles[*it];
            //if(particle -> Status() == -2) {
            //    std::cout << particle -> Position().Pz() << " " << sqrt(radius2) << std::endl;
            //    std::cout << *particle << std::endl;
            //}
            if(particle -> Position().Magnitude2() > radius2 && particle -> Status() != -2) {
                particle -> SetStatus(1);
                it = kickedIdxs.erase(it);
            } else if(particle -> Status() == -2 && particle -> Position().Pz() > sqrt(radius2)) {
                it = kickedIdxs.erase(it);
            } else {
                ++it;
            }
        }

        for(auto particle : particles) {
            if(particle.Position()[0] != particle.Position()[0])
                throw;
        }
    }

    for(auto particle : particles) {
        if(particle.Status() == -1) {
            for(auto p : particles) std::cout << p << std::endl;
            throw std::runtime_error("Cascade has failed. Insufficient max steps.");
        }
    }

    return particles;
}

Particles Cascade::MeanFreePath(
    const Particles& _particles,
    const double& kf,
    const double& _radius2,
    const std::size_t& maxSteps) {

    Particles particles = _particles;
    fermiMomentum = kf;
    radius2 = _radius2;
    auto idx = kickedIdxs[0];
    Particle* kickNuc = &particles[idx];
    
    if (kickedIdxs.size() != 1) {
        std::runtime_error("MeanFreePath: only one particle should be kicked.");
    }
    if (kickNuc -> Status() != -3) {
        std::runtime_error(
            "MeanFreePath: kickNuc must have status -3 "
            "in order to accumulate DistanceTraveled."
            );
    }
    bool hit = false;
    for(std::size_t step = 0; step < maxSteps; ++step) {
        // Are we alrady outside nucleus?
        if (kickNuc->Position().Magnitude2() >= radius2) break;
        AdaptiveStep(particles, distance);
        // Identify nearby particles which might interact
        auto nearby_particles = AllowedInteractions(particles, idx);
        if (nearby_particles.size() == 0) continue;
        // Did we hit?
        auto hitIdx = Interacted(particles, *kickNuc, nearby_particles);
        if (hitIdx == -1) continue;
        Particle* hitNuc = &particles[hitIdx];
        // Did we *really* hit? Finalize momentum, check for Pauli blocking.
        hit = FinalizeMomentum(*kickNuc, *hitNuc);
        // Stop as soon as we hit anything
        if (hit) break;
    }
    if (!hit) kickNuc->SetStatus(1);  // Escaped
    kickedIdxs.clear();  // Clean up
    return particles;
}

void Cascade::AdaptiveStep(const Particles& particles, const double& distance) noexcept {
    double beta = 0;
    for(auto idx : kickedIdxs) {
        if(particles[idx].Beta().Magnitude() > beta)
            beta = particles[idx].Beta().Magnitude();
    }

    timeStep = distance/(beta*Constant::HBARC);
}

bool Cascade::BetweenPlanes(const ThreeVector& position,
                                    const ThreeVector& point1,
                                    const ThreeVector& point2) const noexcept {
    // Get distance between planes
    const ThreeVector dist = point2 - point1;

    // Determine if point is between planes
    return ((position - point1).Dot(dist) >= 0 && (position - point2).Dot(dist) <=0);
}

const ThreeVector Cascade::Project(const ThreeVector& position,
                                                   const ThreeVector& planePt,
                                                   const ThreeVector& planeVec) const noexcept {
    // Project point onto a plane
    const ThreeVector projection = (position - planePt).Dot(planeVec)*planeVec; 
    return position - projection;
}

const InteractionDistances Cascade::AllowedInteractions(Particles& particles,
                                                                        const std::size_t& idx) const noexcept {
    InteractionDistances results;

    // Build planes
    const ThreeVector point1 = particles[idx].Position();
    particles[idx].Propagate(timeStep);
    const ThreeVector point2 = particles[idx].Position();
    auto normedMomentum = particles[idx].Momentum().Vec3().Unit();

    // Build results vector
    for(std::size_t i = 0; i < particles.size(); ++i) {
        // TODO: Should particles propagating be able to interact with 
        //       other propagating particles?
        if (particles[i].Status() < 0) continue;
        // if(i == idx) continue;
        // if(particles[i].InFormationZone()) continue;
        if(!BetweenPlanes(particles[i].Position(), point1, point2)) continue;
        auto projectedPosition = Project(particles[i].Position(), point1, normedMomentum);
        double dist2 = (projectedPosition - point1).Magnitude2();

        results.push_back(std::make_pair(i, dist2));
    }

    // Sort array by distances
    std::sort(results.begin(), results.end(), sortPairSecond);

    return results;
}

const double Cascade::GetXSec(const Particle& particle1,
                                      const Particle& particle2) const {
    return m_interactions -> CrossSection(particle1, particle2);
}

int Cascade::Interacted(const Particles& particles, const Particle& kickedParticle,
                                const InteractionDistances& dists) noexcept {
    for(auto dist : dists) {
        const double xsec = GetXSec(kickedParticle, particles[dist.first]);
        const double prob = exp(-M_PI*dist.second/(xsec/10.));
        if(rng.uniform(0.0, 1.0) < prob) return dist.first;
    }

    return -1;
}

bool Cascade::FinalizeMomentum(Particle& particle1,
                                       Particle& particle2) noexcept {
    // Boost to center of mass frame
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM), p2CM = p2Lab.Boost(-boostCM);

    // Generate outgoing momentum
    bool samePID = particle1.PID() == particle2.PID(); 
    double ecm = (p1CM + p2CM).M();
    const double pcm = particle1.Momentum().Vec3().Magnitude() * Constant::mN / ecm;
    std::array<double, 2> rans;
    rng.generate(rans, 0.0, 1.0);
    ThreeVector momentum = m_interactions -> MakeMomentum(samePID, p1CM.Vec3().Magnitude(), pcm, rans);

    FourVector p1Out = FourVector(momentum[0], momentum[1], momentum[2], p1CM.E()); 
    FourVector p2Out = FourVector(-momentum[0], -momentum[1], -momentum[2], p1CM.E()); 
   
    // Boost back to lab frame
    p1Out = p1Out.Boost(boostCM);
    p2Out = p2Out.Boost(boostCM);

    // Assign momenta to particles
    particle1.SetMomentum(p1Out);
    particle2.SetMomentum(p2Out);

    // Check for Pauli Blocking
    bool hit = !(PauliBlocking(particle1) || PauliBlocking(particle2));

    if(hit) {
        // Assign formation zone
        particle1.SetFormationZone(p1Lab, p1Out);
        particle2.SetFormationZone(p2Lab, p2Out);

        // Hit nucleon is now propagating
        // Users are responsibile for updating the status externally as desired
    } else {
        // Assign momenta to particles
        particle1.SetMomentum(p1Lab);
        particle2.SetMomentum(p2Lab);
    }

    return hit;
}

bool Cascade::PauliBlocking(const Particle& particle) const noexcept {
    if(particle.Status() == -2 && particle.Position().Magnitude2() > radius2) return false;
    return particle.Momentum().Vec3().Magnitude() < fermiMomentum;
}
